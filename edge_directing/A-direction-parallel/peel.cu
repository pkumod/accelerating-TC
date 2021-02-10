#include <omp.h>

#include "buildCSR.cu"
#define thrustSortByKeySize 4000000000;

#define threadNum 256
#define blockNum 10000//100 when under test,10000 when run normal size dataset
using namespace std;
clock_t allStart,allEnd,tStart,tEnd;

void initGPU(uint32_t d) {
	int* warmup = NULL;
    cudaMalloc(&warmup, sizeof(int));
    cudaFree(warmup);
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	if (deviceCount == 0) {
			fprintf(stderr, "error: no devices supporting CUDA.\n");
			exit(EXIT_FAILURE);
	}
	int dev = d;
	cudaSetDevice(dev);

	cudaDeviceProp devProps;
	if (cudaGetDeviceProperties(&devProps, dev) == 0)
	{
			printf("Using device %d:\n", dev);
			printf("%s; global mem: %luB; compute v%d.%d; clock: %d kHz; shared mem: %dB; block threads: %d; SM count: %d\n",
			devProps.name, devProps.totalGlobalMem,
			(int)devProps.major, (int)devProps.minor,
			(int)devProps.clockRate,
			devProps.sharedMemPerBlock, devProps.maxThreadsPerBlock, devProps.multiProcessorCount);
	}
	//	setenv("CUDA_DEVICE_MAX_CONNECTIONS", "32", 1);
	cout<<"GPU selected"<<endl;
	//cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
	return ;
}

void outPutResults(const char *input_file , const char *output_file, 
				   uint32_t *h_flags, uint32_t *h_minorPriority)
{
	fstream fRead(input_file, ios::in|ios::binary);
	ofstream fWrite(output_file, ios::out|ios::binary);
	edge e;
	while (fRead.read((char *)&e, sizeof(edge))) {
		uint32_t src = FIRST(e), dst = SECOND(e);
		int flag_d = (int)h_flags[src] - h_flags[dst];
		int minor_d = (int)h_minorPriority[src] - h_flags[dst];
		if (flag_d > 0 || (flag_d == 0 && minor_d > 0) ||
			(flag_d == 0 && minor_d == 0 && src > dst)) {
			uint32_t tmp = src;
			src = dst;
			dst = tmp;
		}
		e = ((edge)src << 32)|dst;
		fWrite.write((char *)&e, sizeof(edge));
	}
	fRead.close();
	fWrite.close();
	return ;
}		

template <typename T>
struct less_equ {
	T base;
	less_equ(T b) : base(b) {}
	__host__ __device__ bool operator() (const T x) { return x <= base; }
};

template <typename T>
struct my_equal {
	T base;
	my_equal(T b) : base(b) {}
	__host__ __device__ bool operator() (const T x) { return x == base; }
};
__global__ void advance(uint32_t *frontier, uint32_t f_size, uint32_t *flags, 
						uint8_t *visited_flags, int *edgeRow, int *edgeOffset) {
	uint32_t warp_idx = (threadIdx.x + blockDim.x * blockIdx.x)/32;
	uint32_t total_warp = (blockDim.x * gridDim.x)/32;
	uint32_t lane_id = threadIdx.x % 32;
	for (uint32_t i = warp_idx; i < f_size; i += total_warp) {
		uint32_t cur_node = frontier[i];
		uint32_t cur_nbr_size = edgeOffset[cur_node + 1] - edgeOffset[cur_node];
		uint32_t base_off = edgeOffset[cur_node];
		for (uint32_t j = lane_id; j < cur_nbr_size; j += 32) {
			uint32_t cur_nbr = edgeRow[base_off + j];
			if (flags[cur_nbr] == 0xFFFFFFFF)
				visited_flags[cur_nbr] = 1;
		}
	}
	return ;
}
__global__ void back_nbr_count(uint32_t *frontier, uint32_t f_size, uint8_t *last_frontier_flags, 
							   int *edgeRow, int *edgeOffset, uint32_t *runtimeDegree) {
	uint32_t warp_idx = (threadIdx.x + blockDim.x * blockIdx.x)/32;
	uint32_t total_warp = (blockDim.x * gridDim.x)/32;
	uint32_t lane_id = threadIdx.x % 32;
	uint32_t sum = 0;
	for (uint32_t i = warp_idx; i < f_size; i += total_warp) {
		uint32_t cur_node = frontier[i];
		uint32_t cur_nbr_size = edgeOffset[cur_node + 1] - edgeOffset[cur_node];
		uint32_t base_off = edgeOffset[cur_node];
		sum = 0;
		for (uint32_t j = lane_id; j < cur_nbr_size; j += 32) {
			uint32_t cur_nbr = edgeRow[base_off + j];
			if (last_frontier_flags[cur_nbr] != 0)
				sum ++;
		}
		sum += __shfl_down(sum,16);
		sum += __shfl_down(sum, 8);
		sum += __shfl_down(sum, 4);
		sum += __shfl_down(sum, 2);
		sum += __shfl_down(sum, 1);
		if (lane_id == 0) runtimeDegree[cur_node] -= sum;
	}
	return ;
}

__global__ void template_copy(uint32_t *frontier, uint32_t  f_size,
						 uint32_t *runtimeDegree, uint32_t *minorPriority) {
	uint32_t idx = threadIdx.x + blockDim.x * blockIdx.x;
	uint32_t total_thread = blockDim.x * gridDim.x;
	for (uint32_t i = idx; i < f_size; i += total_thread) {
		uint32_t cur_node = frontier[i];
		minorPriority[cur_node] = runtimeDegree[cur_node];
	}
	return ;
}
int main(int argc, const char * argv[]) {
	if (argc != 3) {
		cout << "Usage: ./Peel inputFileName outFileName" << endl;
		return 0;
	}
	initGPU(3);

	tStart = clock();
	unsigned int nodeNum;
	unsigned int edgeNum;
	if (!preProcess(argv[1], edgeNum, nodeNum)) {
		cout << "preprocess failed!" << endl;
		return 0;
	}
	tEnd = clock();
	cout << "preWork cost " << (double)1000*(tEnd-tStart)/CLOCKS_PER_SEC << " ms." << endl;
	cout << "the node num is " << nodeNum << ", and the edgeNum is " << edgeNum << endl;
	/*move csr to GPU**********************************************************************/
	int *d_edgeOffset;
	int *d_edgeRow;
	cudaMalloc(&d_edgeOffset, sizeof(unsigned int)*(nodeNum+2));
	cudaMalloc(&d_edgeRow, sizeof(int)*(1+edgeNum));
	cudaMemcpy((void *)d_edgeOffset, (void *)edgeOffset, sizeof(unsigned int)*(nodeNum+2), cudaMemcpyHostToDevice);
	cudaMemcpy((void *)d_edgeRow, (void *)edgeRow, sizeof(int)*(1+edgeNum), cudaMemcpyHostToDevice);
	cout << "move csr to GPU done!" << endl;
	/***********************************************************************/
	float deg_thre = (float)edgeNum/2/nodeNum;
	//for Debug
	//float deg_thre = 2.2;
	uint32_t total_visited_vertex = 0;
	uint32_t *d_flags;
	uint8_t *d_visited_flags;
	cudaMalloc((void **)&d_visited_flags, sizeof(uint8_t)*nodeNum);
	cudaMalloc((void **)&d_flags, sizeof(uint32_t)*nodeNum);//use this as major priority
	cudaMemset(d_flags, -1, sizeof(uint32_t)*nodeNum);
	uint32_t *d_runtimeDegree;
	uint32_t *d_minorPriority;
	cudaMalloc((void **)&d_minorPriority, sizeof(uint32_t)*nodeNum);//TODO:initialization
	cudaMalloc((void **)&d_runtimeDegree, sizeof(uint32_t)*nodeNum);
	cudaMemcpy(d_runtimeDegree, degreeRecord, sizeof(uint32_t)*nodeNum, cudaMemcpyHostToDevice);
	cudaMemcpy(d_minorPriority, degreeRecord, sizeof(uint32_t)*nodeNum, cudaMemcpyHostToDevice);
	uint32_t *frontier;
	cudaMalloc((void **)&frontier, sizeof(uint32_t)*nodeNum);
	double ss = omp_get_wtime();
	//collect first frontier
	typedef thrust::device_ptr<uint32_t> tptr;
	thrust::device_vector<uint32_t> seq(nodeNum);
	thrust::device_vector<uint32_t> values(nodeNum);
	thrust::device_vector<uint32_t> max_values(nodeNum, 0xFFFFFFFF);
	thrust::device_vector<uint8_t> last_frontier_flag(nodeNum, 0);
	thrust::sequence(thrust::device, seq.begin(), seq.end());
	tptr thr_runtimeDegree(d_runtimeDegree);
	tptr thr_frontier(frontier);
	tptr thr_minorPriority(d_minorPriority);
	tptr index_end = thrust::copy_if(seq.begin(), seq.end(),thr_runtimeDegree,
								   thr_frontier, less_equ<uint32_t>((uint32_t)deg_thre));
	uint32_t f_size = index_end - thr_frontier;
	total_visited_vertex += f_size;
	//set flag to 1 of elements in frontier
	thrust::fill(thrust::device, values.begin(), values.end(), 1);
	thrust::scatter(thrust::device, values.begin(), values.begin()+f_size,
									thr_frontier, 
									thrust::device_pointer_cast(d_flags));
	std::cout << "the frontier size is " << f_size << std::endl;
	//collect runtime degree as edge direction minor priority
	template_copy<<<blockNum, threadNum>>>(frontier, f_size,
								   		   d_runtimeDegree,d_minorPriority);
	cudaDeviceSynchronize();
	//set the d_runtimeDegree to -1 for elements in frontier for later frontier collection
	thrust::scatter(thrust::device, max_values.begin(), max_values.begin() + f_size,
								    thr_frontier,
									thr_runtimeDegree);
	//iterate over increasing deg_thre
	uint32_t peel_depth = 1;
	//cout << "Debug the degree threshold is " << deg_thre << endl;
	//for Debug
	/*uint32_t *h_frontier = new uint32_t [nodeNum];
	cudaMemcpy(h_frontier, frontier, sizeof(uint32_t)*f_size, cudaMemcpyDeviceToHost);
	cout << "frontier: ";
	for (int i = 0; i < f_size; i ++)
		cout << h_frontier[i] << ",";
	cout << endl;
	uint32_t *h_uint32_vec = new uint32_t [nodeNum];
	cudaMemcpy(h_uint32_vec, d_flags, sizeof(uint32_t)*nodeNum, cudaMemcpyDeviceToHost);
	cout << "h_flags after the first frontier" << endl;
	for (int i = 0; i < nodeNum; i++)
		cout << unsigned(h_uint32_vec[i]) << ",";
	cout << endl;*/
	std::cout << "enter the loop here" << std::endl;
	uint32_t itera_times = 1;
	while (total_visited_vertex < nodeNum*9/10 && peel_depth < 0xFF && itera_times < 50) {
		itera_times ++;
		cudaMemset(d_visited_flags, 0, sizeof(uint8_t)*nodeNum);
		//advance from current frontier
		advance<<<blockNum, threadNum>>>(frontier, f_size, d_flags, d_visited_flags, 
										 d_edgeRow, d_edgeOffset);
		cudaDeviceSynchronize();
		//make last frontier flags
		thrust::fill(thrust::device, last_frontier_flag.begin(), last_frontier_flag.end(), 0);
		thrust::scatter(thrust::device, values.begin(), values.begin()+f_size,
										thr_frontier,
										last_frontier_flag.begin());
		//collect reached frontier
		tptr end_pos = thrust::copy_if(seq.begin(),seq.end(),thrust::device_pointer_cast(d_visited_flags),
									   thr_frontier, my_equal<uint8_t>(1));
		//calulate how much the degree declines for vertices in last frontier, write in d_runtimedegree
		f_size = end_pos - thr_frontier;
		back_nbr_count<<<blockNum, threadNum>>>(frontier, f_size, 
												thrust::raw_pointer_cast(last_frontier_flag.data()), 
												d_edgeRow, d_edgeOffset,d_runtimeDegree);
		cudaDeviceSynchronize();
		//collect next peeling frontier
		end_pos = thrust::copy_if(thrust::device, seq.begin(), seq.end(),thr_runtimeDegree,
								  thr_frontier, less_equ<uint32_t>((uint32_t)deg_thre*peel_depth));
		f_size = end_pos - thr_frontier;
		while (f_size == 0 && peel_depth < 0xFF && total_visited_vertex < 9*nodeNum/10) {
			peel_depth ++;
			end_pos = thrust::copy_if(thrust::device, seq.begin(), seq.end(),thr_runtimeDegree,
								     thr_frontier, less_equ<uint32_t>((uint32_t)deg_thre*peel_depth));
			f_size = end_pos - thr_frontier;
			cout << "the f_size equals 0, so we increase the peel_depth as " << peel_depth << endl;
		}
		cout << "the f_size now is " << f_size << endl;
		total_visited_vertex += f_size;
		//set d_flags	
		thrust::fill(thrust::device, values.begin(), values.end(), itera_times);
		thrust::scatter(thrust::device, values.begin(), values.begin()+f_size,
										thr_frontier, 
										thrust::device_pointer_cast(d_flags));
		//collect runtime degree as edge direction minor priority
		template_copy<<<blockNum,threadNum>>>(frontier,f_size,
								   	   		  d_runtimeDegree,d_minorPriority);
		cudaDeviceSynchronize();
		//set the d_runtimeDegree to -1 for elements in frontier for later frontier collection
		thrust::scatter(thrust::device, max_values.begin(), max_values.begin() + f_size,
								    	thr_frontier,
										thr_runtimeDegree);
		if (f_size <= nodeNum/500) peel_depth ++;
		cout << "itera is " << itera_times << ", peel_depth is " << peel_depth 
			  << ", f_size is " << f_size << ", the total_visited_vertex is "
				<< total_visited_vertex << endl;
		//for Debug
		/*cudaMemcpy(h_frontier, frontier, sizeof(uint32_t)*f_size, cudaMemcpyDeviceToHost);
		cout << "frontier: ";
		for (int i = 0; i < f_size; i ++)
			cout << h_frontier[i] << ",";
		cout << endl;
		uint32_t *h_flag = new uint32_t [nodeNum];
		int *h_runtimeDegree = new int [nodeNum];
		cudaMemcpy(h_flag, d_flags, sizeof(uint32_t)*nodeNum, cudaMemcpyDeviceToHost);
		cudaMemcpy(h_runtimeDegree, d_runtimeDegree, sizeof(uint32_t)*nodeNum,cudaMemcpyDeviceToHost);
		cout << "here is the h_flag: " << endl;
		for (int i = 0; i < nodeNum; i ++) 
			cout << unsigned(h_flag[i]) << " ";
		cout << endl << "here is the h_runtimeDegree: " << endl;
		for (int i = 0; i < nodeNum; i ++)
			cout << h_runtimeDegree[i] << " ";
		cout << endl << endl;
		delete [] h_flag;
		delete [] h_runtimeDegree;*/
	}
	//for remained vertices, we simple use "runtimeDegree" as "minorPriority"

	double ee = omp_get_wtime();

	uint32_t *h_flags = new uint32_t [nodeNum];
	uint32_t *h_minorPriority = new uint32_t [nodeNum];
	cudaMemcpy(h_flags, d_flags, sizeof(uint32_t)*nodeNum, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_minorPriority, d_minorPriority, sizeof(uint32_t)*nodeNum, cudaMemcpyDeviceToHost);
	outPutResults(argv[1], argv[2], h_flags, h_minorPriority);
	delete [] h_flags;
	delete [] h_minorPriority;

	cudaFree(d_minorPriority);
	cudaFree(d_visited_flags);
	cudaFree(d_flags);
	cudaFree(d_runtimeDegree);
	cudaFree(frontier);

	delete [] edgeOffset;
	delete [] edgeRow;
	delete [] degreeRecord;
	cudaFree(d_edgeRow);
	cudaFree(d_edgeOffset);
	//cout << "kernel use " << (ee - startKernel) << " s." << endl;
	cout << "The edge diretion uses " << (ee -ss) << " s." << endl;
	return 0;
}
