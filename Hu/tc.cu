#include <omp.h>

#include "preTC.cu"
#define thrustSortByKeySize 4000000000;

#define threadsPerBlockInTC 256
#define blocksPerKernelInTC 30000//100 when under test,10000 when run normal size dataset
#define shareMemorySizeInBlock 3000
#define totalBitsInSharedMemory 96000
#define nodesATimeInBlock 1000
#define maxVNodeATimeInBlock 50000
//#define workPerThreadB 30;
//#define workPerThreadM 20;
//#define maxThreadsPerEdge 256;
using namespace std;
__constant__ unsigned int * c_offset;
__constant__ int * c_row; 
__constant__ int * c_adjLen;
__constant__ int * c_blockStartNodeOffset;
__constant__ long int * c_sum;
clock_t allStart,allEnd,tStart,tEnd;
__global__ void getNodesWorkLoad(int startPos, int threadNum, long int *d_nodeWorkLoad){
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (idx >= threadNum)
		return;
	int src = startPos + idx;
	int *srcList = c_row + c_offset[src];
	int srcListLen = c_adjLen[src];
	long int totalLength = 0;
	for (int i = 0; i < srcListLen; i ++) {
		totalLength += c_adjLen[srcList[i]];
	}
	d_nodeWorkLoad[idx] = totalLength*(unsigned)log2((double)srcListLen);
	return ;
}
unsigned int binarySearchValue(long int *array, long int value, unsigned int arrayLength, int direction) {
	long int s = 0, e = arrayLength-1;
	long int rightPos;
	bool find = false;
	long int mid = (s + e)/2;
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
			if (e == 0)
				break;
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

__global__ void triangleCountKernel() {
	int uid = c_blockStartNodeOffset[blockIdx.x];
	int uidThre = c_blockStartNodeOffset[blockIdx.x+1];
	if (uid == uidThre) return ;
	unsigned vpos = c_offset[uid];
	int vid = c_row[vpos];
	int wpos = threadIdx.x;
	//int *srcList = c_row + c_offset[uid];
	//int *dstList = c_row + c_offset[c_row[vpos]];
	long int sum = 0;
	__shared__ int cachedMaxUid;
	__shared__ bool cacheWorked;
	__shared__ int sharedVid[shareMemorySizeInBlock];
	__shared__ int lastCachedMaxUid;
	while (1) {
		__syncthreads();
		if (threadIdx.x == 0) {
			if (uid == c_blockStartNodeOffset[blockIdx.x]) {
				cachedMaxUid = uid;
			}
			lastCachedMaxUid = cachedMaxUid;
			cachedMaxUid --;
			cacheWorked = true;
			int cachedVid = 0;
			int cachedUid = 0;
			while (1) {
				cachedMaxUid ++;
				cachedUid ++;
				cachedVid += c_adjLen[cachedMaxUid];
				if (cachedVid >= shareMemorySizeInBlock || cachedMaxUid >= uidThre) {
					break;
				}
			}
			if (cachedUid == 1) {
				cacheWorked = false;
				cachedMaxUid++;
				cachedMaxUid = (uidThre > cachedMaxUid)? cachedMaxUid: uidThre;
			} else {
				int len = cachedVid - c_adjLen[cachedMaxUid];
				memcpy(sharedVid,c_row+c_offset[lastCachedMaxUid],sizeof(int)*len);
			}
		}
		/*if (cacheWorked) {
			for (int i = 0; i < (shareMemorySizeInBlock+blockDim.x-1)/blockDim.x; i++) {
				int pos = i*threadsPerBlockInTC + threadIdx.x;
				if (pos >= shareMemorySizeInBlock)
					continue;
				sharedVid[pos] = 0;
			}
		}*/
		/*__syncthreads();
		if (cacheWorked) {
			int _uid = lastCachedMaxUid;
			int _vpos = threadIdx.x + c_offset[_uid];
			while (1) {
				while (_vpos >= c_offset[_uid+1]) {
					_uid++;
					if (_uid >= cachedMaxUid)
						break;
					//printf("_vpos is %d, c_adjLen[_uid] is %d, _uid is %d\n",_vpos,c_adjLen[_uid],_uid); 
				}
				if (_uid >= cachedMaxUid)
					break;
				int hashedValue = (murMurHash(&_uid,4)+c_row[_vpos])%totalBitsInSharedMemory;
				//printf("_uid is %d, hashed is %u, c_row[vpos] is %d, hashedValue is %d\n",_uid, murMurHash(&_uid,4), c_row[vpos],hashedValue);
				atomicOr(validBits+(hashedValue>>5),1<<(hashedValue&31));
				_vpos += blockDim.x;
			}
		}*/
		//TODO : only hash vid
		__syncthreads();
				
		while (vpos >= c_offset[uid+1])
			uid ++;
		while(1) {
			while (wpos >= c_adjLen[vid]) {
					wpos -= c_adjLen[vid];
					vpos ++;
					vid = c_row[vpos];
					while (vpos >= c_offset[uid+1]) {
							uid ++;
					}
			}
			if (uid >= cachedMaxUid) {			
					break;
			}
			int * dstList = c_row + c_offset[vid];
			int targetValue = dstList[wpos];
			if (!cacheWorked) {
					int * srcList = c_row + c_offset[uid];
					int s = 0, e = c_adjLen[uid];
					int mid = (s+e)>>1;
					while (s+1 < e) {
							if (srcList[mid] <= targetValue) {
									s = mid;
							} else {
									e = mid;
							}
							mid = (s+e)>>1;
					}
					if (srcList[s] == targetValue)
							sum++;
			} else {
					int adjListOff = c_offset[uid] - c_offset[lastCachedMaxUid];
					int s = 0, e = c_adjLen[uid];
					int mid = (s+e)>>1;
					while (s+1 < e) {
							if (sharedVid[adjListOff+mid] <= targetValue) {
									s = mid;
							} else {
									e = mid;
							}
							mid = (s+e)>>1;
					}
					if (sharedVid[adjListOff+s] == targetValue)
							sum++;
			}
			wpos += threadsPerBlockInTC;
		}
		if (cachedMaxUid >= uidThre)
			break;
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
	if (argc != 3) {
		cout << "Usage: ./TC -f inputFileName" << endl;
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
	int dev = 2;
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
	if (!preProcess(argv[2], edgeNum, nodeNum)) {
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
	cout << "move csr to GPU done!" << endl;
	/*get work load of each node**********************************************************************/
	long int *nodeWorkLoad = new long int [nodeNum];
	//unsigned int nodeWorkLoad[nodeNum];
	int maxThreadsPerBlock = 1024;
	int maxBlocksPerGrid = 64000;
	int maxThreadsPerKernel = maxBlocksPerGrid*maxThreadsPerBlock;
	long int *d_nodeWorkLoad;
	cudaMalloc(&d_nodeWorkLoad,sizeof(long int)*maxThreadsPerKernel);
	for (int i = 0; i < (nodeNum+maxThreadsPerKernel-1)/maxThreadsPerKernel; i ++) {
		int curThread = maxThreadsPerKernel;
		int remainedNodes = nodeNum - i*maxThreadsPerKernel;
		curThread = (remainedNodes > curThread) ? curThread : remainedNodes;
		cudaMemset(d_nodeWorkLoad,0,sizeof(long)*curThread);
		getNodesWorkLoad<<<(curThread+maxThreadsPerBlock-1)/maxThreadsPerBlock,maxThreadsPerBlock>>>(i*maxThreadsPerKernel,curThread,d_nodeWorkLoad);
		cudaMemcpy(nodeWorkLoad+i*maxThreadsPerKernel,d_nodeWorkLoad,sizeof(long int)*curThread,cudaMemcpyDeviceToHost);
	}
	/*
	cudaFree(d_nodeWorkLoad);
	cudaMalloc(&d_nodeWorkLoad,sizeof(unsigned int)*nodeNum);
	cudaMemcpy((void *)d_nodeWorkLoad,(void *)nodeWorkLoad,sizeof(unsigned int)*nodeNum,cudaMemcpyHostToDevice);
	thrust::inclusive_scan(nodeWorkLoad,nodeWorkLoad+nodeNum,nodeWorkLoad);
	cudaMemcpy((void *)nodeWorkLoad,(void *)d_nodeWorkLoad,sizeof(unsigned int)*nodeNum,cudaMemcpyDeviceToHost);
	cudaFree(d_nodeWorkLoad);
	*/
	for (int i = 1; i < nodeNum; i++)
		nodeWorkLoad[i] += nodeWorkLoad[i-1];
	//for (int i = 0; i < nodeNum; i++)
	//	cout << "for node " << i << ", the workload is " << nodeWorkLoad[i] << endl;
	cout << "get work load of each node done!" << endl;
	/*get block start node offset**********************************************************************/
	int threadsPerKernelInTC = threadsPerBlockInTC*blocksPerKernelInTC;
	int maxWarpPerGrid = threadsPerKernelInTC/32;
	long int workLoadStep = (nodeWorkLoad[nodeNum-1] + blocksPerKernelInTC - 1)/blocksPerKernelInTC;
	int *blockStartNodeOffset = new int [blocksPerKernelInTC+1];
	blockStartNodeOffset[0] = 0;
#pragma omp parallel for
	for (int i = 1; i < blocksPerKernelInTC; i ++) {
		blockStartNodeOffset[i] = (int)binarySearchValue(nodeWorkLoad,i*workLoadStep+1,nodeNum,-1);
	}
	blockStartNodeOffset[blocksPerKernelInTC] = nodeNum;
	//for (int i = 0; i <= maxWarpPerGrid; i ++) {
	//	cout << warpStartNodeOffset[i] << "    " << nodeWorkLoad[warpStartNodeOffset[i]]<< endl;
	//}
	delete [] nodeWorkLoad;
	cout << "get warp start node offset done " << endl;
	/*launch kernel to get result**********************************************************************/
	int *d_blockStartNodeOffset;
	cudaMalloc(&d_blockStartNodeOffset,sizeof(int)*(blocksPerKernelInTC+1));
	cudaMemcpy((void *)d_blockStartNodeOffset,(void *)blockStartNodeOffset,sizeof(int)*(blocksPerKernelInTC+1),cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_blockStartNodeOffset,&d_blockStartNodeOffset,sizeof(int *));
	long int *d_sum;
	cudaMalloc(&d_sum,sizeof(long int)*maxWarpPerGrid);
	cudaMemset(d_sum,0,sizeof(long int)*maxWarpPerGrid);
	cudaMemcpyToSymbol(c_sum,&d_sum,sizeof(long int *));
	double startKernel = omp_get_wtime();
	triangleCountKernel<<<blocksPerKernelInTC,threadsPerBlockInTC>>>();
	cudaDeviceSynchronize();
	//triangleCount = thrust::reduce((thrust::device_ptr<long>)d_sum,(thrust::device_ptr<long>)(d_sum + maxWarpPerGrid));
	long *sum = new long [maxWarpPerGrid];
	cudaMemcpy((void *)sum,(void *)d_sum,sizeof(long)*maxWarpPerGrid,cudaMemcpyDeviceToHost);
	
	triangleCount = thrust::reduce(sum,sum + maxWarpPerGrid);
	delete [] sum;
	cudaFree(d_sum);
	/***********************************************************************/
	/***********************************************************************/
	delete [] blockStartNodeOffset;
	delete [] edgeOffset;
	delete [] edgeRow;
	delete [] adjLength;
	cudaFree(d_blockStartNodeOffset);
	cudaFree(d_edgeRow);
	cudaFree(d_edgeOffset);
	cudaFree(d_adjLength);
	double ee = omp_get_wtime();
	cout << "There are " << triangleCount << " triangles in the input graph." << endl;
	cout << "kernel use " << (ee - startKernel) << " s." << endl;
	cout << "Total use time " << (ee -ss) << " s." << endl;
	return 0;
}
