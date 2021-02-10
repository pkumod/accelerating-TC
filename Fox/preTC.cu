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
#include<thrust/execution_policy.h>
using namespace std;
#define thrustSortBlockSize 4000000000
#define bucketNum 10
struct edge{
    unsigned int src;
    unsigned int dst;
};

struct cmpStruc{
    __device__ bool operator () (const edge &a, const edge &b){
        return (a.src < b.src) || (a.src == b.src && a.dst < b.dst) ;
    }
}cmp;
class edgeVector{
		public:
				unsigned int capcity;
				unsigned int esize;
				edge *Edges;
				edgeVector(){esize = 0; capcity = 0;}
				void init(unsigned int s) { Edges = new edge [s]; capcity = s; return ;}
				void addEdge(int _src,int _dst){
						if(esize >= capcity) {
								capcity *= 2;
								edge* tmpEdges = new edge [capcity];
								memcpy(tmpEdges,Edges,sizeof(edge)*esize);
								delete [] Edges;
								Edges = tmpEdges;
						}
						Edges[esize].src = _src;
						Edges[esize].dst = _dst;
						esize ++;
				}
				void clear() {delete [] Edges; return ;}
};
int *edgeOffset;
int *edgeRow;
int *adjLength;
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
    Edges = new edge [edgeNum];

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

    //pre work
	
	unsigned int maxNodeID = 0;
	//#pragma omp parallel for
    for(unsigned int i = 0; i < edgeNum; i ++)
    {
    	if (Edges[i].src < Edges[i].dst) {
    		int tmpValue = Edges[i].dst;
    		Edges[i].dst = Edges[i].src;
    		Edges[i].src = tmpValue;
    	}
		if (Edges[i].src > maxNodeID)
			maxNodeID = Edges[i].src;
    }
	cout << "end rearrange dst and src" << endl;
	unsigned nodeNum = maxNodeID + 1;
    //sort edges
	//************sort edges && get nodeNum********
 	edgeVector * edgeBucket = new edgeVector [bucketNum];
	for (int i = 0; i < bucketNum; i ++) 
		edgeBucket[i].init(edgeNum/bucketNum);
	unsigned bucketStep = (nodeNum + bucketNum - 1)/bucketNum; 
	for (int i = 0; i < edgeNum; i ++)
	{
		int bucketID = Edges[i].src/bucketStep;
		edgeBucket[bucketID].addEdge(Edges[i].src, Edges[i].dst);
	}
	cout << "end pust edges in bucket" << endl;
	unsigned int *bucketEdgeOffset = new unsigned int [bucketNum];
	bucketEdgeOffset[0] = 0;
	for (int i = 0; i < bucketNum-1; i ++) {
		unsigned int bucketSize = edgeBucket[i].esize;
		if (bucketSize > thrustSortBlockSize/sizeof(edge)) {
			cout << "bucket " << i << "size is " << bucketSize << ", it's too large!" << endl;
			return false;
		}
		bucketEdgeOffset[i+1] = bucketEdgeOffset[i] + bucketSize;
	}
	for (int i = 0; i < bucketNum; i++) {
		thrust::device_vector<edge> D (edgeBucket[i].Edges, edgeBucket[i].Edges+edgeBucket[i].esize);
		thrust::sort(D.begin(),D.begin()+edgeBucket[i].esize,cmp);
		thrust::copy(D.begin(),D.begin()+edgeBucket[i].esize,edgeBucket[i].Edges);
	}
	cout << "end sort edges in GPU " << endl;
	for(int i = 0; i < bucketNum; i ++) {
		memcpy(Edges+bucketEdgeOffset[i],edgeBucket[i].Edges,sizeof(edge)*edgeBucket[i].esize);
	}
	cout << "end copy result to Edges" << endl;
	delete [] bucketEdgeOffset;
	for (int i = 0; i < bucketNum; i ++)
		edgeBucket[i].clear();
	delete [] edgeBucket;
	//************end sort edges && get nodeNum********

    //unsigned int nodeNum = Edges[edgeNum-1].src + 1;
    edgeOffset = new int [nodeNum+1];
    edgeRow = new int [edgeNum];
    adjLength = new int[nodeNum];
    unsigned int nodePos = 0;
    unsigned int edgePos = 0;
    edgeOffset[0] = 0;
    edge * edgePtr;
    int formerSrc = -1,formerDst = -1;
	start_ = clock();
    for (unsigned int i = 0; i < edgeNum; i++)
    {
	    edgePtr = Edges + i;
//	    cout << "start put cur edge in csr" << endl;
	    if (edgePtr->src == edgePtr->dst) {
		    formerSrc = edgePtr->src;
		    formerDst = edgePtr->dst;
		    continue;
	    }
//	    cout << "1 " << endl;
	    if ((i > 0) && (edgePtr->src == formerSrc)) {
		    //TODO  find a more efficienty way
		    if(edgePtr->dst == formerDst){
			    continue;
		    }
		    edgeRow[edgePos++] = edgePtr->dst;
		    formerDst = edgePtr->dst;
		    continue;
	    }
//	    cout << "2 " << endl;
	    int curSrc = edgePtr->src;
	    for (unsigned j = nodePos + 1; j <= curSrc; j++) {
		    edgeOffset[j] = edgePos;
		    adjLength[j-1] = edgeOffset[j]-edgeOffset[j-1];	
	    }
	    nodePos = curSrc;
	    edgeRow[edgePos++] = edgePtr->dst;
	    formerSrc = edgePtr->src;
	    formerDst = edgePtr->dst;
//	    cout << " end an edge in a loop " << endl;
    }
	end_ = clock();
	cout << "merge and make csr use " << (double)1000*(end_-start_)/CLOCKS_PER_SEC << " ms." << endl;
    edgeOffset[nodeNum] = edgePos;
    adjLength[nodeNum-1] = edgeOffset[nodeNum] - edgeOffset[nodeNum-1];
    cout << "csr built, edgeNum is "<< edgePos<< ", the node num is " << nodeNum << ", origin egde num is " << edgeNum << endl;
    //TODO remove empty node in edgeOffset
    _edgeNum = edgeOffset[nodeNum];
    _nodeNum = nodeNum;
    delete [] Edges;
	cout << "rebuild Edges now " << endl;
	double start_omp = omp_get_wtime();
	Edges = new edge[_edgeNum];

	#pragma omp parallel for
	for (int i = 0; i < _nodeNum; i++) {
		int *curList = edgeRow + edgeOffset[i];
		for (int j = 0; j < adjLength[i]; j++) {
				Edges[edgeOffset[i]+j].src = i;
				Edges[edgeOffset[i]+j].dst = curList[j];
		}
	}
	double end_omp = omp_get_wtime();
	cout << "rebuild use " << (end_omp-start_omp) << " s."<< endl;
	cout << "rebuild done" << endl;
    return true;
}
