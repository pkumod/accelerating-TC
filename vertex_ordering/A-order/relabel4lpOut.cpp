#include<iostream>
#include<string.h>
#include<fstream>
#include<stdlib.h>
#include<vector>
#include<algorithm>
#include<random>
using namespace std;
struct Edge{
	int src;
	int dst;
};
int main(int argc, const char * argv[]){
	fstream fRead;
	fRead.open(argv[1],ios::in);
	int maxNodeId = 0;
	Edge curEdge;
	vector<pair<int,int> > allEdges;
	int src,dst;
	while(fRead >> src >> dst) {
		maxNodeId = (src > maxNodeId) ? src : maxNodeId;
		maxNodeId = (dst > maxNodeId) ? dst : maxNodeId;
		allEdges.push_back(make_pair(src,dst));
	}
	fRead.close();

	int *degree = new int [maxNodeId+1];
	memset(degree,0,sizeof(int)*(maxNodeId+1));
	for (int i = 0; i < allEdges.size(); i ++) {
		degree[allEdges[i].first]++;
	}

	int maxRunTimeDegree = 0;
	for (int i = 0; i <= maxNodeId; i++) {
		if (degree[i] > maxRunTimeDegree)
			maxRunTimeDegree = degree[i];
	}
	cout << "the max runtime degree is " << maxRunTimeDegree << endl;
	vector<int> *bucket = new vector<int> [maxRunTimeDegree+1];
	for (int i = 0; i <= maxNodeId; i ++) {
		bucket[degree[i]].push_back(i);
	}
	vector<pair<int,int> > id_rtDegree;
	for (int i = 0; i <= maxRunTimeDegree; i ++){
		for (int j = 0; j < bucket[i].size(); j++) {
			id_rtDegree.push_back(make_pair(bucket[i][j],i));
		}
	}
	
	int *goodLookUpTable = new int [maxNodeId+1];
	memset(goodLookUpTable,0,sizeof(int)*(maxNodeId+1));
	//TODO
	for (int i = 0; i < id_rtDegree.size(); i ++) {
		goodLookUpTable[id_rtDegree[i].first] = id_rtDegree.size()-1-i;
	}
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	int mixRatio = 40;
	int * badDegreeDistri = new int [maxNodeId+1];
	int startPos = 0;
	memset(badDegreeDistri,0,sizeof(int)*(maxNodeId+1));
	//cout << 1 << endl;
	//*******************************************method 1 start
	/*int cutNum = (maxNodeId+mixRatio)/mixRatio;
	for (int i = 0; i < cutNum; i++) {
		badDegreeDistri[startPos] = id_rtDegree[maxNodeId-i].first;
		startPos += mixRatio;
	}
	//cout << 2 << endl;
	startPos = 1;
//	for (int i = (maxNodeId+1)/2; i <= maxNodeId; i ++) {
	for (int i = 0; i <= maxNodeId-cutNum; ) {
		for (int k = 0; k < mixRatio-1; k ++) {
			badDegreeDistri[startPos++] = id_rtDegree[i++].first;
			if (i > maxNodeId-cutNum)
				break;
		}
		startPos++;
	}*/
	//******************************************method 1 done
	//******************************************method 2 start
	
	reverse(id_rtDegree.begin(),id_rtDegree.end());
	bool reverse = false;
	int chunkSize = (maxNodeId+1)/mixRatio;
	for (int i = 0; i < mixRatio; i ++) {
		if(!reverse){
			int startWritePos = i;
			int startReadPos = i*chunkSize;
			for (int j = 0; j < chunkSize; j ++) {
				badDegreeDistri[startWritePos] = id_rtDegree[startReadPos].first;
				startWritePos += mixRatio;
				startReadPos ++;
			}
		} else {
			int startWritePos = i;
			int startReadPos = (i+1)*chunkSize-1;
			for (int j = 0; j < chunkSize; j ++) {
				badDegreeDistri[startWritePos] = id_rtDegree[startReadPos].first;
				startWritePos += mixRatio;
				startReadPos --;
			}
		}
		reverse = !reverse;
	}
	for (int i = 0; i < chunkSize; i ++) {
		random_shuffle(badDegreeDistri+i*mixRatio, badDegreeDistri+(i+1)*mixRatio);
	}
	for (int i = (1+maxNodeId)-(1+maxNodeId)%mixRatio; i <= maxNodeId; i ++) {
		badDegreeDistri[i] = id_rtDegree[i].first;
	}


	int * badLookUpTable = new int [maxNodeId+1];
	memset(badLookUpTable,0,sizeof(int)*(maxNodeId+1));
	for (int i = 0; i <= maxNodeId; i++) {
		badLookUpTable[badDegreeDistri[i]] = i;
	}

	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//	for (int i = 0; i < maxNodeId; i ++) 
//		cout << runTimeDegree[id_rtDegree[i].first] << endl;
	ofstream fWrite_normal("normal",ios::out|ios::binary);
	ofstream fWrite_bad("bad",ios::out|ios::binary);
	ofstream fWrite_good("good",ios::out|ios::binary);
	for (int i = 0; i < allEdges.size(); i ++) {
		int src = allEdges[i].first;
		int dst = allEdges[i].second;
		curEdge.src = src;
		curEdge.dst = dst;
		fWrite_normal.write((char *)&curEdge,sizeof(Edge));
		curEdge.src = goodLookUpTable[src];
		curEdge.dst = goodLookUpTable[dst];
		fWrite_good.write((char*)&curEdge,sizeof(Edge));
		curEdge.src = badLookUpTable[src];
		curEdge.dst = badLookUpTable[dst];
		fWrite_bad.write((char *)&curEdge,sizeof(Edge));
	}
	fWrite_normal.close();
	fWrite_good.close();
	fWrite_bad.close();
	delete [] degree;
	allEdges.clear();
	for (int i = 0; i <= maxRunTimeDegree; i ++)
		bucket[i].clear();
	delete [] bucket;
	delete [] goodLookUpTable;
	delete [] badDegreeDistri;
	delete [] badLookUpTable;
	id_rtDegree.clear();
	return 0;
}
