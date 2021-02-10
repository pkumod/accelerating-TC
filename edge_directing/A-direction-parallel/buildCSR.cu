#include <iostream>
#include <algorithm>
#include <vector>
#include <iterator>
#include <fstream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <cmath>
#include<sys/stat.h>
#include<ctime>

#include <cuda_runtime.h>

#include<thrust/reduce.h>
#include<cuda_runtime.h>
#include<thrust/sort.h>
#include<thrust/device_ptr.h>
#include<thrust/device_vector.h>
#include<thrust/host_vector.h>
#include<thrust/copy.h>
#include<thrust/fill.h>
#include<thrust/scatter.h>
#include<thrust/execution_policy.h>
#include<thrust/scan.h>
#include<thrust/iterator/counting_iterator.h>
#include<thrust/functional.h>
#include<iterator>
using namespace std;
#define thrustSortBlockSize 4000000000
#define bucketNum 20

typedef uint64_t edge;
#define FIRST(x) (uint32_t(x>>32u))
#define SECOND(x) (uint32_t(x))

struct is_self_loop : public thrust::unary_function<uint64_t, bool> {
	__host__ __device__ bool operator() (uint64_t x) {
		return FIRST(x) == SECOND(x);
	}
};
struct cmpStruc{
    __device__ bool operator () (const edge &a, const edge &b){
        return (FIRST(a) < FIRST(b)) || (FIRST(a) == FIRST(b) && SECOND(a) < SECOND(b)) ;
    }
}cmp;

class edgeVector{
	public:
		unsigned int capcity;
		unsigned int esize;
		edge *Edges;
		edgeVector(){esize = 0; capcity = 0;}
		void init(unsigned int s) { Edges = new edge [s]; capcity = s; return ;}
		void addEdge(edge * E){
			if(esize >= capcity) {
				capcity *= 2;
				edge* tmpEdges = new edge [capcity];
				memcpy(tmpEdges,Edges,sizeof(edge)*esize);
				delete [] Edges;
				Edges = tmpEdges;
			}
			memcpy(Edges+esize,E,sizeof(edge));
			esize ++;
		}
		void clear() {delete [] Edges; return ;}
};
unsigned int *edgeOffset;
int *edgeRow;
int *degreeRecord;
edge *Edges;
clock_t start_, end_;
bool preProcess(const char *fileName, unsigned int  &_edgeNum, unsigned &_nodeNum)
{
    //get file size
	
    ifstream fin1(fileName,ios::in|ios::binary);
    fin1.seekg(0,ios::end);
    streampos Size = fin1.tellg();
    fin1.close();
    long int size = Size;
    cout << "the size of input file is " << size << " Byte. " << endl;
    unsigned int edgeNum = size/(sizeof(int)*2);
    Edges = new edge [edgeNum*2];

    //read data
    ifstream fin(fileName, std::ios::binary);
    if (fin.bad()) {
		cout << "File not fould!" << endl;
		return false;
	}
    cout << "start read data... ..." << endl;
    fin.read((char *)Edges,sizeof(edge)*edgeNum);
    fin.close();
    cout << "end read data" << endl;

    //find node number
	uint32_t nodeNum = 0;
#pragma omp parallel for reduction(max : nodeNum)
    for (uint64_t i = 0; i < edgeNum; ++i) {
        edge e = Edges[i];
        nodeNum = max(nodeNum, FIRST(e));
        nodeNum = max(nodeNum, SECOND(e));
    }
    nodeNum++;

	//cal degrees
	degreeRecord = new int[nodeNum];
	memset(degreeRecord,0,sizeof(int)*nodeNum);
	//#############################################
	//duplicate edges
	#pragma omp parallel for
	for (unsigned int i = 0; i < edgeNum; i ++) {
		unsigned int src = FIRST(Edges[i]);
		unsigned int dst = SECOND(Edges[i]);
		Edges[i+edgeNum] = ((edge)dst << 32) | src;
	}
	edgeNum *= 2;	
	
	cout << "end rearrange dst and src" << endl;
    //sort edges
 	edgeVector * edgeBucket = new edgeVector [bucketNum];
	for (int i = 0; i < bucketNum; i ++) 
		edgeBucket[i].init(edgeNum/bucketNum);
	unsigned bucketStep = (nodeNum + bucketNum - 1)/bucketNum; 
	for (unsigned int i = 0; i < edgeNum; i ++)
	{
		int bucketID = FIRST(Edges[i])/bucketStep;
		edgeBucket[bucketID].addEdge(Edges+i);
	}
	cout << "end pust edges in bucket" << endl;
	for (int i = 0; i < bucketNum-1; i ++) {
		unsigned int bucketSize = edgeBucket[i].esize;
		if (bucketSize > thrustSortBlockSize/sizeof(edge)) {
			cout << "bucket " << i << "size is " << bucketSize << ", it's too large!" << endl;
			return false;
		}
	}
	edgeNum = 0;
	for (int i = 0; i < bucketNum; i++) {
		uint32_t effective_len = edgeBucket[i].esize;
		thrust::device_vector<edge> D (edgeBucket[i].Edges, edgeBucket[i].Edges + edgeBucket[i].esize);
		thrust::sort (D.begin(), D.begin() + edgeBucket[i].esize, cmp);
		effective_len = thrust::remove_if(D.begin(), D.begin() + edgeBucket[i].esize, is_self_loop()) - D.begin();
        effective_len = thrust::unique(D.begin(), D.begin() + effective_len) - D.begin();
		thrust::copy(D.begin(), D.begin()+ effective_len, Edges + edgeNum);
		edgeNum += effective_len;
	}
	cout << "end sort edges in GPU and copy result to Edges" << endl;
	
	for (int i = 0; i < bucketNum; i ++)
		edgeBucket[i].clear();
	delete [] edgeBucket;
	
    edgeOffset = new unsigned int [nodeNum+1];
    edgeOffset[0] = 0;
    edgeRow = new int [edgeNum];

	start_ = clock();

	int node_pos = FIRST(Edges[0]);
	for (int i = 0; i <= node_pos; i ++)
		edgeOffset[i] = 0;
	edgeRow[0] = SECOND(Edges[0]);
	for (unsigned int i = 1; i < edgeNum; i++)
	{
		edge cur_edge = Edges[i];
		edgeRow[i] = SECOND(cur_edge);
		int src = FIRST(cur_edge);
		if (src > node_pos) {
			for (int j = node_pos + 1; j <= src; j ++)
				edgeOffset[j] = i;
			node_pos = src;
		}
    }
	for (int i = node_pos + 1; i <= nodeNum; i ++) {
			edgeOffset[i] = edgeNum;
	}

	end_ = clock();

	for (int i = 0; i < nodeNum; i ++)
		degreeRecord[i] = edgeOffset[i+1] - edgeOffset[i];

	cout << "merge and make csr use " << (double)1000*(end_-start_)/CLOCKS_PER_SEC << " ms." << endl;
    cout << "csr built, edgeNum is "<< edgeNum << ", the node num is " << nodeNum << ", origin egde num is " << edgeNum << endl;
    //TODO remove empty node in edgeOffset
    _edgeNum = edgeNum;
    _nodeNum = nodeNum;
    delete [] Edges;
	//CSR validation check
	/*for (int i = 0 ; i < nodeNum; i ++) {
		cout << "the neighbor list of vertex " << i << ", size is :" << degreeRecord[i] << endl;
		for (int j = edgeOffset[i]; j < edgeOffset[i+1]; j ++)
			cout << edgeRow[j] << " ";
		cout << endl;
	}*/
    return true;
}
