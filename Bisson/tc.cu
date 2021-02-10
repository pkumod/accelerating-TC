#include <omp.h>
#include <cassert>
#include "preTC.cu"

#define shareMemorySizeInBlock 200
#define hIndex 2048

using namespace std;
__constant__ unsigned int * c_offset;
__constant__ int * c_row; 
__constant__ int * c_adjLen;
__constant__ long int * c_sum;
__constant__ int * c_bitmap;
__constant__ int * c_nonZeroRow;
__device__ int nextNode;
clock_t allStart,allEnd,tStart,tEnd;

__global__ void triangleCountKernel(unsigned int totalNodeNum, int nonZeroSize) {
	long int sum = 0;
	int curRowNum = blockIdx.x;
	int lane_id = threadIdx.x%32;
	__shared__ int sh_bitMap[shareMemorySizeInBlock];
	unsigned int intSizePerBitmap = (totalNodeNum+31)/32;
	int * myBitmap = c_bitmap + blockIdx.x*intSizePerBitmap;
	while(1) {
		//__syncthreads();
		int privateRowNum = (curRowNum < nonZeroSize)? c_nonZeroRow[curRowNum] : totalNodeNum;
		//int privateRowNum = curRowNum;
		if (privateRowNum >= totalNodeNum) {
			break;
		}
		//if (c_offset[privateRowNum+1] == c_offset[privateRowNum])
		//	continue;
		int * curNodeNbr = c_row + c_offset[privateRowNum];
		unsigned int curNodeNbrLength = c_offset[privateRowNum+1] - c_offset[privateRowNum];
		if (curNodeNbrLength > 32) {
				if (threadIdx.x == 0) {
						memset(myBitmap, 0,sizeof(int)*intSizePerBitmap);
						memset(sh_bitMap,0,sizeof(int)*shareMemorySizeInBlock);
				}
				__threadfence();
				for (int i = 0; i < (curNodeNbrLength+blockDim.x-1)/blockDim.x; i ++) {
						int curIndex = i*blockDim.x + threadIdx.x;
						int curNbr;
						if (curIndex < curNodeNbrLength) {
								curNbr = curNodeNbr[curIndex];
								atomicOr(myBitmap+(curNbr/32),1<<(31-curNbr%32));
								atomicOr(sh_bitMap+(curNbr/hIndex/32),1<<(31-(curNbr/hIndex)%32));
						}
						__syncthreads();
						if (curIndex < curNodeNbrLength) {
								int * twoHoopNbr = c_row + c_offset[curNbr];
								unsigned int  twoHoopNbrLength = c_offset[curNbr+1] - c_offset[curNbr];
								for (int j = 0; j < twoHoopNbrLength; j ++) {
										int curValue = twoHoopNbr[j];
										if (  ((sh_bitMap[curValue/hIndex/32]>>(31-(curValue/hIndex)%32))&1)
														&&((myBitmap[curValue/32]>>(31-curValue%32))&1) ) {
												sum++;
										}
								}
						}
				}
		} else {
			if (threadIdx.x == 0)
				memcpy(sh_bitMap, curNodeNbr, sizeof(int)*curNodeNbrLength);
			__threadfence();
			for (int i = lane_id; i < curNodeNbrLength; i += 32) {
					int curNbr = curNodeNbr[i];
					int *twoHoopNbr = c_row + c_offset[curNbr];
					int twoHoopNbrLength = c_offset[curNbr+1] - c_offset[curNbr];
					for (int j = 0; j < twoHoopNbrLength; j ++) {
						int targetValue = twoHoopNbr[j];
						int s = 0, e = curNodeNbrLength, mid;
						while (s < e) {
							mid = (s + e)/2;
							if (sh_bitMap[mid] > targetValue)
								e = mid;
							else if (sh_bitMap[mid] < targetValue)
								s = mid + 1;
							else{
								sum ++;
								break;
							}
						}	
					}
			}
		}
		curRowNum += gridDim.x;//atomicAdd(&nextNode, 1);
		//__syncthreads();
		//if (privateRowNum != curRowNum) 
		//	printf("private is %d, curRowNum is %d, block %d, index %d\n",privateRowNum,curRowNum,blockIdx.x,threadIdx.x);
	}		

	sum += __shfl_down(sum,16);
	sum += __shfl_down(sum,8);
	sum += __shfl_down(sum,4);
	sum += __shfl_down(sum,2);
	sum += __shfl_down(sum,1);
	int idx = blockDim.x*blockIdx.x + threadIdx.x;
	if (threadIdx.x%32 == 0) {
			c_sum[idx>>5] = sum;
	}
	return ;
}
int main(int argc, const char * argv[]) {
	/**********************************prework of the algorithm: read data & make CSR*********************************/
	if (argc != 4) {
		cout << "Usage: ./TC -f inputFileName chooseindex" << endl;
		return 0;
	}
	
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
	//cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

	tStart = clock();
	double ss = omp_get_wtime();
	unsigned int nodeNum;
	unsigned int edgeNum;
	int nonZeroSize;
	if (!preProcess(argv[2], edgeNum, nodeNum, nonZeroSize, atoi(argv[3]))) {
		cout << "preprocess failed!" << endl;
		return 0;
	}
	tEnd = clock();
	cout << "preWork cost " << (double)1000*(tEnd-tStart)/CLOCKS_PER_SEC << " ms." << endl;
	cout << "the node num is " << nodeNum << ", and the edgeNum is " << edgeNum << endl;
	
	long int triangleCount = 0;
	/*move csr to GPU**********************************************************************/
	int *d_adjLength;
	cudaMalloc(&d_adjLength,sizeof(int)*(1+nodeNum));
	cudaMemcpy((void *)d_adjLength,(void *)adjLength,sizeof(int)*(1+nodeNum),cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_adjLen,&d_adjLength,sizeof(int *));
	int *d_edgeOffset;
	int *d_edgeRow;
	cudaMalloc(&d_edgeOffset, sizeof(unsigned int)*(nodeNum+2));
	cudaMalloc(&d_edgeRow, sizeof(int)*(1+edgeNum));
	cudaMemcpy((void *)d_edgeOffset, (void *)edgeOffset, sizeof(unsigned int)*(nodeNum+2), cudaMemcpyHostToDevice);
	cudaMemcpy((void *)d_edgeRow, (void *)edgeRow, sizeof(int)*(1+edgeNum), cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(c_offset,&d_edgeOffset,sizeof(unsigned int *));
	cudaMemcpyToSymbol(c_row,&d_edgeRow,sizeof(int *));

	int *d_nonZeroRow;
	cudaMalloc(&d_nonZeroRow, sizeof(int)*nonZeroSize);
	cudaMemcpy((void*)d_nonZeroRow,(void*)nonZeroRow,sizeof(int)*nonZeroSize,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_nonZeroRow,&d_nonZeroRow,sizeof(int*));

	int h_nextNode = 0;
	cudaMemcpyToSymbol(nextNode,&h_nextNode,sizeof(int));
	cout << "move csr to GPU done!" << endl;

	int bitPerInt = sizeof(int)*8;
	unsigned intSizePerBitmap = (nodeNum+bitPerInt-1)/bitPerInt;
	int blockSize = 32;
	int blockNum = 30*2048/blockSize;
	
	if (nodeNum > hIndex*shareMemorySizeInBlock*32) {
		cout << "ERROR! the nodeNum is too large: " << nodeNum << endl;
		return 0;
	}
	if (blockNum*intSizePerBitmap*sizeof(int)/1024 > 8*1024*1024) {
		cout << "RUN OUT OF GLOBAL MEMORY!!" << endl;
		return 0;
	}
	int *d_bitmaps;
	cudaMalloc(&d_bitmaps,sizeof(int)*intSizePerBitmap*blockNum);
	cudaMemcpyToSymbol(c_bitmap, &d_bitmaps, sizeof(int *));
	/*launch kernel to get result**********************************************************************/
	long int *d_sum;
	unsigned maxWarpPerGrid = blockNum*blockSize/32;
	cudaMalloc(&d_sum,sizeof(long int)*maxWarpPerGrid);
	cudaMemset(d_sum,0,sizeof(long int)*maxWarpPerGrid);
	cudaMemcpyToSymbol(c_sum,&d_sum,sizeof(long int *));
	double startKernel = omp_get_wtime();
	triangleCountKernel<<<blockNum,blockSize>>>(nodeNum,nonZeroSize);
	cudaDeviceSynchronize();
	//triangleCount = thrust::reduce((thrust::device_ptr<long>)d_sum,(thrust::device_ptr<long>)(d_sum + maxWarpPerGrid));
	long *sum = new long [maxWarpPerGrid];
	cudaMemcpy((void *)sum,(void *)d_sum,sizeof(long)*maxWarpPerGrid,cudaMemcpyDeviceToHost);
	
	triangleCount = thrust::reduce(sum,sum + maxWarpPerGrid);
	delete [] sum;
	cudaFree(d_sum);
	/***********************************************************************/
	//for debug 
	/*int testData[] = {34523,34114,34115,34116,34117,2051};
	for (int j = 0; j < 6; j ++) {
		cout << "the nbr of " << testData[j] << " is " << endl;
		for (int i = 0; i < edgeOffset[testData[j]+1]-edgeOffset[testData[j]]; i ++) {
			cout << edgeRow[edgeOffset[testData[j]]+i] << endl;
		}
	}*/

	delete [] edgeOffset;
	delete [] edgeRow;
	delete [] adjLength;
	delete [] nonZeroRow;
	cudaFree(d_edgeRow);
	cudaFree(d_edgeOffset);
	cudaFree(d_adjLength);
	cudaFree(d_bitmaps);
	cudaFree(d_nonZeroRow);
	double ee = omp_get_wtime();
	cout << "There are " << triangleCount << " triangles in the input graph." << endl;
	cout << "kernel use " << (ee - startKernel) << " s." << endl;
	cout << "Total use time " << (ee -ss) << " s." << endl;
	return 0;
}
