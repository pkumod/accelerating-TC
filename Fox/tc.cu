#include <omp.h>

#include "preTC.cu"
#define thrustSortByKeySize 4000000000;
//#define workPerThreadB 30;
//#define workPerThreadM 20;
//#define maxThreadsPerEdge 256;
using namespace std;
__constant__ int * c_offset;
__constant__ int * c_row; 
__constant__ int * c_adjLen;
__constant__ unsigned int c_edgeSize;
__constant__ edge *c_Edges;
__constant__ int c_threadsPerEdge;
__constant__ long *c_sums;
//but largeAdjNodeThrehold/blocksize is possible loop times for node with neighborhood near largeAdjNodeThrehold
clock_t allStart,allEnd,tStart,tEnd;
__global__ void getEdgeWorkLoad(int threadNum, int *d_edgeWorkLoad, edge *d_edgeForWL){
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx >= threadNum)
		return;
	int srcLen = c_adjLen[d_edgeForWL[idx].src];
	int dstLen = c_adjLen[d_edgeForWL[idx].dst];
	//printf("src is %d, dst is %d, the srcLen is %d, the dstLen is %d\n",d_edgeForWL[idx].src,d_edgeForWL[idx].dst,srcLen,dstLen);
	int large = (srcLen > dstLen)?srcLen:dstLen;
	int bWork = (srcLen + dstLen - large)*log2((double)large+2);
	//if (k == 0)
	//	printf("!!!! idx is %d, src degree is %d, dst degree is %d\n",idx,srcLen,dstLen);
	//else 
	//	printf("idx is %d, src degree is %d, dst degree is %d, K = %d\n",idx,srcLen,dstLen,k);
	d_edgeWorkLoad[idx] = bWork;
	return ;
}
unsigned int binarySearchValue(int *array, int value, unsigned int arrayLength, int direction) {
	unsigned int s = 0, e = arrayLength-1;
	unsigned int rightPos;
	bool find = false;
	unsigned int mid = (s + e)/2;
	while (s <= e)
	{
		if(array[mid] == value) {
			rightPos = mid;
			find = true;
			break;
		}
		else if(array[mid] < value) {
			s = mid + 1;
		}
		else {
			e = mid - 1;
		}
		mid = (s + e)/2;
	}
	if (!find) {
		return s;
	}
	long int tmpValue = rightPos + direction;
	while(tmpValue >= 0 && tmpValue < arrayLength && array[tmpValue] == value) {
			rightPos += direction;
			tmpValue = rightPos + direction;
	}
	return rightPos;
}
__global__ void binSearchKernel(){
	long idx = blockDim.x*blockIdx.x + threadIdx.x;
	long sum = 0;
	__shared__ long sh_sum[32];
	//sh_sum[threadIdx.x/32] = 0;
	int edgeID = idx/c_threadsPerEdge;
	int inEdgeID = idx%c_threadsPerEdge;
	
	while (edgeID < c_edgeSize) { 
		int src = c_Edges[edgeID].src;
		int dst = c_Edges[edgeID].dst;	
		if (c_adjLen[src] < c_adjLen[dst]) {
			int tmp = src;
			src = dst;
			dst = tmp;
		}
		int srcAdjListLen = c_adjLen[src];
		int dstAdjListLen = c_adjLen[dst];
		int *srcAdjList = c_row + c_offset[src];
		int *dstAdjList = c_row + c_offset[dst];
		for (int i = 0; i < (dstAdjListLen+c_threadsPerEdge-1)/c_threadsPerEdge; i++) {
			int dstListIdx = i*c_threadsPerEdge + inEdgeID;
			if (dstListIdx >= dstAdjListLen) 
				continue;
			int targetValue = dstAdjList[dstListIdx];
			int s = 0,e = srcAdjListLen-1;	
			int mid = (s + e)/2;
			while (s <= e)
			{
					if(srcAdjList[mid] == targetValue) {
							sum++;
							break;
					}
					else if(srcAdjList[mid] < targetValue) {
							s = mid + 1;
					}
					else {
							e = mid - 1;
					}
					mid = (s + e)/2;
			}
		}
		idx += blockDim.x*gridDim.x;
		edgeID = idx/c_threadsPerEdge;
		inEdgeID = idx%c_threadsPerEdge;
	}
	int tIdx = threadIdx.x;
	for(int i = 16; i >= 1; i /= 2)
		sum += __shfl_down(sum,i);
	if(tIdx%32 == 0)
		sh_sum[tIdx/32] = sum;
	__syncthreads();
	if (tIdx == 0) {
		sum = 0;
		for (int i = 0; i < 32; i++) {
			sum += sh_sum[i];
		}
		c_sums[blockIdx.x] = sum;
	}
	return;

}
double totalReorderTime;
void rearrangeEdge(edge *Edges,
			unsigned int binaryStartPos,
				unsigned curSize,
					 int edgesInABlock,
					int *adjLength ) {
	
	int *srcLen = new int [curSize];
#pragma omp parallel for
	for (int i = 0; i < curSize; i ++) {
		edge *curEdge = Edges + binaryStartPos + i;
		int length1 = adjLength[curEdge->src];
		int length2 = adjLength[curEdge->dst];
		srcLen[i] = length1 > length2 ? length1 : length2;
	}
	edge *curEdges = new edge[curSize];
	memcpy(curEdges,Edges+binaryStartPos,sizeof(edge)*curSize);
	thrust::sort_by_key(srcLen,srcLen+curSize,curEdges);
//	memcpy(Edges+binaryStartPos,curEdges,sizeof(edge)*curSize);

	int extraEdgesNum = curSize%edgesInABlock;
	memcpy(Edges+binaryStartPos+curSize-extraEdgesNum,curEdges,sizeof(edge)*extraEdgesNum);
	int groupSize = 32;
	double start_t = omp_get_wtime();
/*	unsigned i,desPos;
//#pragma omp parallel for
	for (i = curSize-1, desPos = 0; i >= extraEdgesNum; i-= groupSize, desPos += edgesInABlock) {
		if (desPos >= curSize-extraEdgesNum) {
			desPos = desPos%edgesInABlock + groupSize;
		}
//#pragma omp parallel for
		//for (int j = 0; j < groupSize; j++) {
		//	Edges[binaryStartPos+desPos+j].src = curEdges[i-j].src;
		//	Edges[binaryStartPos+desPos+j].dst = curEdges[i-j].dst;
			//i --;
		//}
		memcpy(Edges+binaryStartPos+desPos,curEdges+i-groupSize-1,sizeof(edge)*groupSize);
		//i -= groupSize;
		//desPos += edgesInABlock;
	}*/
	int blockNum = curSize/edgesInABlock;
	int groupPerBlock = edgesInABlock/groupSize;
#pragma omp parallel for
	for (int i = 0; i < curSize/groupSize; i ++) {
		int firstIndex,secIndex;
		firstIndex = i%blockNum;
		secIndex = i/blockNum;
		memcpy(Edges+binaryStartPos+(firstIndex*groupPerBlock+secIndex)*groupSize,
			   curEdges+i*groupSize,
			   sizeof(edge)*groupSize);
	}
	double end_t = omp_get_wtime();
	totalReorderTime += (end_t - start_t);
	cout << "cur reorder time is " << end_t - start_t << endl;
	delete [] curEdges;
	delete [] srcLen;
	return ;
}
int main(int argc, const char * argv[]) {
	/**********************************prework of the algorithm: read data & make CSR*********************************/
	if (argc != 4) {
		cout << "Usage: ./TC -f inputFileName" << endl;
		return 0;
	}
	totalReorderTime = 0.0;
	int* warmup = NULL;
    cudaMalloc(&warmup, sizeof(int));
    cudaFree(warmup);
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	if (deviceCount == 0) {
			fprintf(stderr, "error: no devices supporting CUDA.\n");
			exit(EXIT_FAILURE);
	}
	int dev = 1;
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


	tStart = clock();
	double ss = omp_get_wtime();
	unsigned int nodeNum;
	unsigned int edgeNum;
	if (!preProcess(argv[2], edgeNum, nodeNum)) {
		cout << "preprocess failed!" << endl;
		return 0;
	}
	tEnd = clock();
	cout << "preWork cost " << (double)1000*(tEnd-tStart)/CLOCKS_PER_SEC << " ms." << endl;
	cout << "the node num is " << nodeNum << ", and the edgeNum is " << edgeNum << endl;
	cudaEvent_t CstartTime, CendTime;
	cudaEventCreate(&CstartTime);
	cudaEventCreate(&CendTime);
	//sort & copy data to GPU
	long int triangleCount = 0;
	//TODO move this to preProcess
	/***************************************sort edges according to work load*********************************/
	int *edgeWorkLoad = new int [edgeNum];
	//for (int i = 0; i < nodeNum; i ++)
	//		cout << "node[" << i << "]'s degree is " << adjLength[i] << endl;
	int *d_adjLength;
	cudaMalloc(&d_adjLength,sizeof(int)*nodeNum);
	cudaMemcpy((void *)d_adjLength,(void *)adjLength,sizeof(int)*nodeNum,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_adjLen,&d_adjLength,sizeof(int *));
	cout << "calculate the work load of every edge" << endl;		
	//calculate the work load of every edge 
	int maxBlockNumWL = 64000;
	int maxBlockSizeWL = 1024;
	int maxTotalThreadWL = maxBlockNumWL*maxBlockSizeWL;
	int * d_edgeWorkLoad;
	edge * d_edgesForWL;
	//TODO : add if
	cudaMalloc(&d_edgesForWL,sizeof(edge)*maxTotalThreadWL);
	cudaMalloc(&d_edgeWorkLoad,sizeof(int)*maxTotalThreadWL);
	for (unsigned int i = 0; i < (edgeNum+maxTotalThreadWL-1)/maxTotalThreadWL; i++) {
		int leftEdge = edgeNum - i*maxTotalThreadWL;
		cout << "leftEdge is " << leftEdge << endl;
		int curThread = (leftEdge > maxTotalThreadWL)? maxTotalThreadWL : leftEdge;
		cout << "curThread is " << curThread << endl;
		cudaMemcpy((void *)d_edgesForWL,(void *)(Edges + i*maxTotalThreadWL),sizeof(edge)*curThread,cudaMemcpyHostToDevice);
		int curGridSize = (curThread+maxBlockSizeWL-1)/maxBlockSizeWL;
		cout << "curGridSize is " << curGridSize << endl;
		getEdgeWorkLoad<<<curGridSize,maxBlockSizeWL>>>(curThread,d_edgeWorkLoad,d_edgesForWL);
		cudaDeviceSynchronize();
		cudaMemcpy((void *)(edgeWorkLoad + i*maxTotalThreadWL),(void *)d_edgeWorkLoad,sizeof(int)*curThread,cudaMemcpyDeviceToHost);
	}
	cudaFree(d_edgeWorkLoad);
	cudaFree(d_edgesForWL);
	//for (int i = 0; i < edgeNum; i ++)
	//	cout << "workLoad is " << edgeWorkLoad[i] << ",and the src and dst degree are " << adjLength[Edges[i].src] << " " << adjLength[Edges[i].dst] << endl;
	cout << "sort edges" << endl;
	//sort edges
	unsigned int thrustSortByKeyThre = thrustSortByKeySize;
	thrustSortByKeyThre/=(sizeof(int)*2);
	int totalPass = (edgeNum + thrustSortByKeyThre - 1)/thrustSortByKeyThre;
	edge *d_Edges;
	int *d_edgeWL;
	//TODO :add if
	cout << "totalPass is " << totalPass << endl;
	cudaMalloc(&d_Edges,sizeof(edge)*thrustSortByKeyThre);
	cudaMalloc(&d_edgeWL,sizeof(int)*thrustSortByKeyThre);
	for(int i = 0; i < totalPass; i ++) {
		long int leftEdge = edgeNum - i*thrustSortByKeyThre;
		long int curEdge = (leftEdge > thrustSortByKeyThre) ? thrustSortByKeyThre : leftEdge;
		cudaMemcpy((void *)d_Edges,(void *)(Edges+i*thrustSortByKeyThre),sizeof(edge)*curEdge,cudaMemcpyHostToDevice);
		cudaMemcpy((void *)d_edgeWL,(void *)(edgeWorkLoad+i*thrustSortByKeyThre),sizeof(int)*curEdge,cudaMemcpyHostToDevice);
//		thrust::device_ptr<edge> dev_data_ptr(raw_data);
//		thrust::device_ptr<int> dev_keys_ptr(raw_keys);
		thrust::sort_by_key((thrust::device_ptr<int>)d_edgeWL,(thrust::device_ptr<int>)(d_edgeWL+curEdge),(thrust::device_ptr<edge>)d_Edges);
		cudaMemcpy((void *)(Edges+i*thrustSortByKeyThre),(void *)d_Edges,sizeof(edge)*curEdge,cudaMemcpyDeviceToHost);
		cudaMemcpy((void *)(edgeWorkLoad+i*thrustSortByKeyThre),(void *)d_edgeWL,sizeof(int)*curEdge,cudaMemcpyDeviceToHost);
	}
	
	cudaFree(d_Edges);
	cudaFree(d_edgeWL);
//	for (int i = 0; i < edgeNum; i++) {
//		cout << "workLoad is " << edgeWorkLoad[i] << ",and the src and dst degree are " << adjLength[Edges[i].src] << " " << adjLength[Edges[i].dst] << endl;
//	}
	//mv csr to GPU
	cout << "move csr to GPU" << endl;
	int *d_edgeOffset;
	int *d_edgeRow;
	cudaMalloc(&d_edgeOffset, sizeof(int)*(nodeNum+1));
	cudaMalloc(&d_edgeRow, sizeof(int)*edgeNum);
	cudaMemcpy((void *)d_edgeOffset, (void *)edgeOffset, sizeof(int)*(nodeNum+1), cudaMemcpyHostToDevice);
	cudaMemcpy((void *)d_edgeRow, (void *)edgeRow, sizeof(int)*edgeNum, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(c_offset,&d_edgeOffset,sizeof(int *));
	cudaMemcpyToSymbol(c_row,&d_edgeRow,sizeof(int *));
	double totalKernelTimeInMsec = 0;
	double totalTCTimeInMsec = 0;
	for (int p = 0; p < totalPass; p ++) {//for each big block divided by thrust sort
		cout << "Pass " << p << endl;
		unsigned curBlockStartPos = p*thrustSortByKeyThre;
		unsigned curBlockSize = edgeNum - curBlockStartPos;
		curBlockSize = (curBlockSize > thrustSortByKeyThre)?thrustSortByKeyThre:curBlockSize;
		//seperate two parts
		unsigned int binaryStartPos = curBlockStartPos;
	
		int binMaxWork = edgeWorkLoad[binaryStartPos+curBlockSize-1];
		cout << "binMaxWork is " << binMaxWork << endl;
		int workPerThreadB = 8;
		int maxThreadsPerEdge = 32;
		double maxPowerOf8 = log(binMaxWork)/log(workPerThreadB);
		//******************
		int curBinB = maxThreadsPerEdge;
		unsigned int *binStartPos = new unsigned int[curBinB+2];
		binStartPos[1] = 0;
		for(int i = 2; i <= curBinB; i*=2) {
			double index = (double)(maxPowerOf8-1)*(i-1)/(double)curBinB+1;
			int partitionPoint = powl(workPerThreadB,index);
			binStartPos[i] = binarySearchValue(edgeWorkLoad+binaryStartPos,partitionPoint,curBlockSize,-1);
			binStartPos[i/2+1] = binStartPos[i];
		}
		binStartPos[curBinB+1] = curBlockSize;
	//	cout << "binStartPos[" << curBinB << "] is " << binaryLength << endl;
		//#pragma omp parallel for
		for(int i = 1; i <= curBinB; i*=2) {//for each bin
			tStart = clock();
			unsigned int curSize = binStartPos[i+1] - binStartPos[i];
			if (curSize == 0) {
				continue;
			}
			cout<< "start point workload " << edgeWorkLoad[binStartPos[i]] << ", end point workload " << edgeWorkLoad[binStartPos[i+1]-1] << endl;
		//	cout<< "start point " << binStartPos[i] << ", end point" << binStartPos[i+1] << endl;
			int maxBlockNumTC = 20000;
			int maxBlockSizeTC = 1024;
			int curThreadsPerEdge = i;
			if (atoi(argv[3])) {
				rearrangeEdge(Edges+binaryStartPos,binStartPos[i],curSize,maxBlockSizeTC/curThreadsPerEdge,adjLength);
			}
			edge * d_curEdges;
			cudaMalloc(&d_curEdges,sizeof(edge)*curSize);
			long *d_sum;
			cudaMalloc(&d_sum,sizeof(long)*maxBlockNumTC);
			long *h_sum = new long[maxBlockNumTC];
			cout << "curThreadsPerEdge is " << curThreadsPerEdge << endl;
			cout << "curSize is " << curSize << endl;
			cudaMemcpyToSymbol(c_Edges,&d_curEdges,sizeof(edge *));
			cudaMemcpyToSymbol(c_threadsPerEdge,&curThreadsPerEdge,sizeof(int));
			cudaMemcpyToSymbol(c_sums,&d_sum,sizeof(long*));

			cudaMemset(d_sum,0,sizeof(long)*maxBlockNumTC);
			cudaMemcpy((void *)d_curEdges,(void *)(Edges+binaryStartPos+binStartPos[i]),sizeof(edge)*curSize,cudaMemcpyHostToDevice);
			unsigned int curGridSize = 20000;	
			//		cout << "	curGridSize is " << curGridSize << endl;
			cudaMemcpyToSymbol(c_edgeSize,&curSize,sizeof(unsigned int));
			cudaEventRecord(CstartTime,0);
			binSearchKernel<<<curGridSize,maxBlockSizeTC>>>();
			cudaDeviceSynchronize();
			cudaEventRecord(CendTime,0);
			cudaEventSynchronize(CstartTime);
			cudaEventSynchronize(CendTime);
			float time;
			cudaEventElapsedTime(&time,CstartTime,CendTime);
			totalKernelTimeInMsec += time;
			//cout << "one kernel finished! " << endl;
			cudaMemcpy((void *)h_sum,(void *)d_sum,sizeof(long)*curGridSize,cudaMemcpyDeviceToHost);
			triangleCount += thrust::reduce(h_sum,h_sum+curGridSize);

			cudaFree(d_curEdges);
			cudaFree(d_sum);
			delete [] h_sum;
			tEnd = clock();
			cout << " curTriangleNum is " << triangleCount << endl;
			cout << " use time " << time << " ms." << endl << endl;
			totalTCTimeInMsec += (double)1000*(tEnd-tStart)/CLOCKS_PER_SEC;
		}//for each bin
		delete [] binStartPos;

	}//for each big block divided by thrust sort
	/***************************trangle counting on GPU******************************/
	delete [] edgeWorkLoad;
	delete [] edgeOffset;
	delete [] edgeRow;
	delete [] adjLength;
	delete [] Edges;
	cudaFree(d_edgeRow);
	cudaFree(d_edgeOffset);
	cudaFree(d_adjLength);
	cudaEventDestroy(CstartTime);
	cudaEventDestroy(CendTime);
	double ee = omp_get_wtime();
	allEnd = clock();
	cout << "the reorder time is " << totalReorderTime << endl;
	cout << "The total Time of kernel is " << totalKernelTimeInMsec << " ms." << endl;
	cout << "The total Time of TC is " << totalTCTimeInMsec << " ms." << endl;
	cout << "There are " << triangleCount << " triangles in the input graph." << endl;
	cout << "Total use time " << (ee -ss) << " s." << endl;
	return 0;
}
