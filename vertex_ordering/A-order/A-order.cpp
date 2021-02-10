#include<iostream>
#include<string.h>
#include<fstream>
#include<stdlib.h>
#include<vector>
#include<algorithm>
//#include<random>
#include<math.h>
#include<time.h>
#include<cmath>
#include<random>
#include "omp.h"
using namespace std;
#define BUCKETSIZE 40
struct Edge{
	int src;
	int dst;
};

float bandWidthLookUpTable[10] = 
	{0,0,271,310,399,419,458,478,524,539};
float lambdaLookUpTable[10] = 
	{0,0,145.51,121.98,53.26,23.92,16,14,12.94,10.24};
float getBandwidth(int length) {
	if (length < 10)
		return bandWidthLookUpTable[length];
	else {
		return 384.1*pow(length,0.1491);
	}
}
float getEqualBandWidth(int length) {
	if (length < 10)
		return lambdaLookUpTable[length]*sqrt(getBandwidth(length))/sqrt((float)length);
	else 
		return 12*sqrt(getBandwidth(length))/sqrt((float)length);
}
	
struct nodeBucket{
	int id;
	float memorySuperiority;
	std::vector<int> *nodes;
};

struct lessCmp{
	bool operator () (const nodeBucket & nb1, const nodeBucket & nb2) {
		return nb1.memorySuperiority < nb2.memorySuperiority || 
				(nb1.memorySuperiority == nb2.memorySuperiority && nb1.id < nb2.id);
	}
};
struct greatCmp{
	bool operator () (const nodeBucket & nb1, const nodeBucket & nb2) {
		return nb1.memorySuperiority > nb2.memorySuperiority || 
				(nb1.memorySuperiority == nb2.memorySuperiority && nb1.id > nb2.id);
	}
};
struct idSort{
	bool operator () ( const nodeBucket &nb1, const nodeBucket &nb2) {
		return nb1.id > nb2.id;
	}
};

int main(int argc, const char * argv[]){
	fstream fRead;
	fRead.open(argv[1],ios::in|ios::binary);
	int maxNodeId = 0;
	Edge curEdge;
	vector<pair<int,int> > allEdges;
	while(fRead.read((char *)&curEdge,sizeof(Edge))) {
		int src = curEdge.src;
		int dst = curEdge.dst;
		maxNodeId = (src > maxNodeId) ? src : maxNodeId;
		maxNodeId = (dst > maxNodeId) ? dst : maxNodeId;
		allEdges.push_back(make_pair(src,dst));
	}
	fRead.close();


	int *degree = new int [maxNodeId+1];
	memset(degree,0,sizeof(int)*(maxNodeId+1));
	for (int i = 0; i < allEdges.size(); i ++) {
		degree[allEdges[i].first]++;
		degree[allEdges[i].second]++;
	}


	int *runTimeDegree = new int [maxNodeId+1];
	memset(runTimeDegree,0,sizeof(int)*(maxNodeId+1));
	for (int i = 0; i < allEdges.size(); i ++) {
		int src = allEdges[i].first;
		int dst = allEdges[i].second;
		if (degree[src] > degree[dst] || 
			(degree[src] == degree[dst] && src < dst)) {
			allEdges[i].first = dst;
			allEdges[i].second = src;
			runTimeDegree[dst]++;
		}
		else
			runTimeDegree[src]++;
	}
	int maxRunTimeDegree = 0;
	for (int i = 0; i <= maxNodeId; i++) {
		if (runTimeDegree[i] > maxRunTimeDegree)
			maxRunTimeDegree = runTimeDegree[i];
	}
	cout << "the max runtime degree is " << maxRunTimeDegree << endl;
	vector<int> *bucket = new vector<int> [maxRunTimeDegree+1];
	for (int i = 0; i <= maxNodeId; i ++) {
		bucket[runTimeDegree[i]].push_back(i);
	}
	vector<pair<int,int> > id_rtDegree;
	unsigned int memoryDominatePos = 0;
	unsigned int twoDegreeNodePos = 0;
	for (int i = 2; i <= maxRunTimeDegree; i ++){
		if (i == 5)
			memoryDominatePos = id_rtDegree.size();
		for (int j = 0; j < bucket[i].size(); j++) {
			id_rtDegree.push_back(make_pair(bucket[i][j],i));
		}
	}

	//cout << 0 << endl;
	cout << "the id_rtDegree size is " << id_rtDegree.size() << endl;
	int totalBucketNum = id_rtDegree.size()/BUCKETSIZE;
	//reverse(id_rtDegree.begin(),id_rtDegree.end());
	//bucket initialization
	vector<nodeBucket> allBucket(totalBucketNum);
	for (int i = 0; i < totalBucketNum; i ++) {
		allBucket[i].id = i;
		allBucket[i].memorySuperiority = 0.0;
		allBucket[i].nodes = new vector<int>;
	}
	make_heap(allBucket.begin(),allBucket.end(),greatCmp());

	float *bandwidth = new float [maxRunTimeDegree+1];
	float *equalBandwidth = new float [maxRunTimeDegree+1];
	for (int i = 0; i < maxRunTimeDegree+1; i ++) {
		bandwidth[i] = getBandwidth(i);
		equalBandwidth[i] = getEqualBandWidth(i);
	}
	int finishedBucket = 0;
	cout << "the effective uid size is " << id_rtDegree.size() << endl;
	double reorderTime = 0.0;
	clock_t start_t = clock();
	for (int i = BUCKETSIZE*totalBucketNum-1; i >= memoryDominatePos; i--) {
			//TODO
		pop_heap(allBucket.begin(),allBucket.end()-finishedBucket,greatCmp());
		struct  nodeBucket &curBucket = *(allBucket.rbegin()+finishedBucket);
		int curNode = id_rtDegree[i].first;
		int curDegree = id_rtDegree[i].second;
		curBucket.nodes->push_back(curNode);
		curBucket.memorySuperiority += (bandwidth[curDegree] - equalBandwidth[curDegree])/BUCKETSIZE;
		if (curBucket.nodes->size() >= BUCKETSIZE) {
			finishedBucket ++;
		} else {
			push_heap(allBucket.begin(),allBucket.end()-finishedBucket,greatCmp());
		}
	}
	clock_t end_t = clock();
	reorderTime += (double)(end_t - start_t)/CLOCKS_PER_SEC;
	cout << "after put all memory-dominated nodes, the finished bucket number is " << finishedBucket << endl; 
//	for (int i = 0; i < totalBucketNum; i++)
//		cout << allBucket[i].memorySuperiority << endl;
	make_heap(allBucket.begin(),allBucket.end()-finishedBucket,lessCmp());
	start_t = clock();
	for (int i = twoDegreeNodePos; i < memoryDominatePos; i ++) {
		pop_heap(allBucket.begin(),allBucket.end()-finishedBucket,lessCmp());
		struct nodeBucket &curBucket = *(allBucket.rbegin()+finishedBucket);
		int curNode = id_rtDegree[i].first;
		int curDegree = id_rtDegree[i].second;
		curBucket.nodes->push_back(curNode);
		curBucket.memorySuperiority -= (equalBandwidth[curDegree]-bandwidth[curDegree])/BUCKETSIZE;
		if (curBucket.nodes->size() >= BUCKETSIZE) {
			finishedBucket ++;
		} else {
			push_heap(allBucket.begin(),allBucket.end()-finishedBucket,lessCmp());
		}
	}
	end_t = clock();
	reorderTime += (double)(end_t-start_t)/CLOCKS_PER_SEC;
	cout << "after put all nodes, the finished bucket number is " << finishedBucket << endl; 
	sort_heap(allBucket.begin(),allBucket.end(),idSort());
	cout << "the total reorder time is " << reorderTime << endl;
	//for (int i = 0; i < totalBucketNum; i++) 
	// random_shuffle(allBucket[i].nodes->begin(),allBucket[i].nodes->end());

	vector<int> badDegreeDistri;
	for (int i = 0; i < totalBucketNum; i ++) {
		struct nodeBucket &curBucket = allBucket[i];
		for (int j = 0; j < curBucket.nodes->size(); j++)
			badDegreeDistri.push_back(curBucket.nodes->at(j));
	}
	int remainedNode = id_rtDegree.size()%BUCKETSIZE;
	for (int i = totalBucketNum*BUCKETSIZE; i < id_rtDegree.size(); i ++)
		badDegreeDistri.push_back(id_rtDegree[i].first);
	for (int i = 0; i < bucket[1].size(); i ++)
		badDegreeDistri.push_back(bucket[1][i]);
	for (int i = 0; i < bucket[0].size(); i ++)
		badDegreeDistri.push_back(bucket[0][i]);

	cout << "useless nodes number is " << bucket[0].size() + bucket[1].size() << endl;
	for (int i = 0; i < totalBucketNum; i ++) {
		if (allBucket[i].nodes->size() != BUCKETSIZE)
			cout << "the bucket id is " << i << ", and size is " << allBucket[i].nodes->size() << endl;
	}
	//cout << 3 << endl;
	int * badLookUpTable = new int [maxNodeId+1];
	memset(badLookUpTable,0,sizeof(int)*(maxNodeId+1));
	for (int i = 0; i <= maxNodeId; i++) {
		badLookUpTable[badDegreeDistri[i]] = i;
	}
	//cout << 4 << endl;
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//	for (int i = 0; i < maxNodeId; i ++) 
//		cout << runTimeDegree[badDegreeDistri[i]] << endl;
//	for (int i = 0; i < totalBucketNum; i++)
//		cout << allBucket[i].memorySuperiority << endl;
	ofstream fWrite_normal("normal",ios::out|ios::binary);
	ofstream fWrite_bad("reordered",ios::out|ios::binary);
	for (int i = 0; i < allEdges.size(); i ++) {
		int src = allEdges[i].first;
		int dst = allEdges[i].second;
		curEdge.src = src;
		curEdge.dst = dst;
		fWrite_normal.write((char *)&curEdge,sizeof(Edge));
		curEdge.src = badLookUpTable[src];
		curEdge.dst = badLookUpTable[dst];
		fWrite_bad.write((char *)&curEdge,sizeof(Edge));
	}
	fWrite_normal.close();
	fWrite_bad.close();
	for (int i = 0; i < totalBucketNum; i ++) {
		allBucket[i].nodes->clear();
		delete  allBucket[i].nodes;
	}
	delete [] runTimeDegree;
	delete [] degree;
	allEdges.clear();
	for (int i = 0; i <= maxRunTimeDegree; i ++)
		bucket[i].clear();
	delete [] bucket;
	delete [] bandwidth;
	delete [] equalBandwidth;
	delete [] badLookUpTable;
	allBucket.clear();
	badDegreeDistri.clear();
	id_rtDegree.clear();
	return 0;
}
