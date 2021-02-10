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
#include<thrust/scan.h>
using namespace std;
#define thrustSortBlockSize 4000000000
#define bucketNum 10
struct edge{
    int src;
    int dst;
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
    //fine node number
	int divideNum = 100;	
	unsigned int *maxNodeIDs = new unsigned int [divideNum];
	memset(maxNodeIDs,0,sizeof(unsigned int)*divideNum);
	#pragma omp parallel for
	for (int d = 0; d < divideNum; d++) {
			unsigned int step = edgeNum/divideNum;
			unsigned int s = d*step;
			unsigned int e = (d+1)*step;
			if (d == divideNum - 1)
				e = edgeNum;
			for(unsigned int i = s; i < e; i ++)
			{
				if (Edges[i].src > maxNodeIDs[d])
					maxNodeIDs[d] = Edges[i].src;
				if (Edges[i].dst > maxNodeIDs[d])
					maxNodeIDs[d] = Edges[i].dst;
			}
	}
	unsigned int maxNodeID = maxNodeIDs[0];
	for (int i = 1; i < divideNum; i ++)
		if (maxNodeIDs[i] > maxNodeID)
			maxNodeID = maxNodeIDs[i];
	cout << "get max nodeid" << endl;
	unsigned nodeNum = maxNodeID + 1;
	delete [] maxNodeIDs;
	//cal degrees
	int * degreeRecord = new int[nodeNum];
	memset(degreeRecord,0,sizeof(int)*nodeNum);
	//#############################################
	for (unsigned int i = 0; i < edgeNum; i++)
	{
		degreeRecord[Edges[i].src]++;
		degreeRecord[Edges[i].dst]++;
	}	

	#pragma omp parallel for
	for (unsigned int i = 0; i < edgeNum; i ++) {
		unsigned int src = Edges[i].src;
		unsigned int dst = Edges[i].dst;
		if (degreeRecord[src] > degreeRecord[dst] || (degreeRecord[src] == degreeRecord[dst] && src < dst)) {
			Edges[i].src = dst;
			Edges[i].dst = src;
		}
	}
		
	int * toBeMini = new int[nodeNum];
	int wipedEdges = 0;
	memset(toBeMini,0,sizeof(int)*nodeNum);
	int totalMinied = 0;
	for (unsigned int i = 0; i < nodeNum; i ++) {
		if (degreeRecord[i] <= 1) {
			totalMinied ++;
		}
		toBeMini[i] = totalMinied;
	}


	#pragma omp parallen for 
	for (unsigned int i = 0; i < edgeNum; i++) {
		unsigned int src = Edges[i].src;
		unsigned int dst = Edges[i].dst;	
		if (degreeRecord[src] <= 1) {
			Edges[i].src = -1;
			wipedEdges ++;
			continue;
		}
		if (degreeRecord[dst] <= 1) {
			Edges[i].dst = -1;
			wipedEdges ++;
			continue;
		}
		if (src > 0) {
			Edges[i].src = src - toBeMini[src-1];
		}
		if (dst > 0) 
			Edges[i].dst = dst - toBeMini[dst-1];
	}
	nodeNum = nodeNum - totalMinied;
	delete []  toBeMini;
	delete [] degreeRecord;	

	
	cout << "end rearrange dst and src" << endl;
	//######################################
	/*#pragma omp parallel for 
	for (unsigned int i = 0; i < edgeNum; i++) {
		unsigned int src = Edges[i].src;
		unsigned int dst = Edges[i].dst;
		if (src < dst) {
			Edges[i].src = dst;
			Edges[i].dst = src;
		}
	}*/
	//#########################################
    //sort edges
	//************sort edges && get nodeNum********
 	edgeVector * edgeBucket = new edgeVector [bucketNum];
	for (int i = 0; i < bucketNum; i ++) 
		edgeBucket[i].init(edgeNum/bucketNum);
	unsigned bucketStep = (nodeNum + bucketNum - 1)/bucketNum; 
	for (int i = 0; i < edgeNum; i ++)
	{
		if (Edges[i].src == -1)
			continue;
		int bucketID = Edges[i].src/bucketStep;
		edgeBucket[bucketID].addEdge(Edges+i);
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
	edgeNum = edgeNum - wipedEdges;//************************************************
    //unsigned int nodeNum = Edges[edgeNum-1].src + 1;
    edgeOffset = new unsigned int [nodeNum+2];
    edgeOffset[0] = 0;
    edgeRow = new int [edgeNum+1];
    adjLength = new int[nodeNum+1];
	memset(adjLength,0,sizeof(int)*(nodeNum+1));
	unsigned int nodePos = 0;
	unsigned int edgePos = 0;
	edge * edgePtr;
	int formerSrc = -1,formerDst = -1;
	start_ = clock();
//	for (int i = 0; i < edgeNum; i++)
//		printf("%d   %d\n",Edges[i].src,Edges[i].dst);
	for (unsigned int i = 0; i < edgeNum; i++)
	{
			edgePtr = Edges + i;
			if (edgePtr->src == -1 || edgePtr->dst == -1)
				continue;
			if (edgePtr->src == edgePtr->dst) {
					formerSrc = edgePtr->src;
					formerDst = edgePtr->dst;
					int curSrc = edgePtr->src;
					for (unsigned j = nodePos + 1; j <= curSrc; j++) {
							edgeOffset[j] = edgePos;
							adjLength[j-1] = edgeOffset[j]-edgeOffset[j-1];	
					}
					nodePos = curSrc;
					continue;
			}
			if ((i > 0) && (edgePtr->src == formerSrc)) {
					//TODO  find a more efficienty way
					if(edgePtr->dst == formerDst){
							continue;
					}
					edgeRow[edgePos++] = edgePtr->dst;
					formerDst = edgePtr->dst;
					continue;
			}	
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
	for (unsigned i = nodePos + 1; i < nodeNum; i ++) {
			edgeOffset[i] = edgePos;
			adjLength[i-1] = edgeOffset[i] - edgeOffset[i-1];
	}

	end_ = clock();
	cout << "merge and make csr use " << (double)1000*(end_-start_)/CLOCKS_PER_SEC << " ms." << endl;
    edgeOffset[nodeNum] = edgePos;
	edgeOffset[nodeNum+1] = edgePos + 1;
    adjLength[nodeNum-1] = edgeOffset[nodeNum] - edgeOffset[nodeNum-1];
	adjLength[nodeNum] = 1024;
	edgeRow[edgePos] = nodeNum;
    cout << "csr built, edgeNum is "<< edgePos<< ", the node num is " << nodeNum << ", origin egde num is " << edgeNum << endl;
    //TODO remove empty node in edgeOffset
	int maxDegreeStored = 0;
	for (int i = 0; i < nodeNum; i ++)
		if (adjLength[i] > maxDegreeStored)
			maxDegreeStored = adjLength[i];
	cout << "The max stored degree is " << maxDegreeStored  << endl;
    _edgeNum = edgeOffset[nodeNum];
    _nodeNum = nodeNum;
    delete [] Edges;
    return true;
}
