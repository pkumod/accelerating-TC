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
	for (int i = 0; i <= maxRunTimeDegree; i ++){
		for (int j = 0; j < bucket[i].size(); j++) {
			id_rtDegree.push_back(make_pair(bucket[i][j],i));
		}
	}
	//cout << -1 << endl;
	int *goodLookUpTable = new int [maxNodeId+1];
	memset(goodLookUpTable,0,sizeof(int)*(maxNodeId+1));
	//TODO
	for (int i = 0; i < id_rtDegree.size(); i ++) {
		goodLookUpTable[id_rtDegree[i].first] = id_rtDegree.size()-1-i;
	}
	//cout << 0 << endl;
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	int mixRatio = 40;
	int * badDegreeDistri = new int [maxNodeId+1];
	for (int i = 0; i <= maxNodeId; i ++)
		badDegreeDistri[i] = i;

	random_shuffle(badDegreeDistri, badDegreeDistri+(maxNodeId+1));
	
	//********************************************method 2 done
	//cout << 3 << endl;
	int * badLookUpTable = new int [maxNodeId+1];
	memset(badLookUpTable,0,sizeof(int)*(maxNodeId+1));
	for (int i = 0; i <= maxNodeId; i++) {
		badLookUpTable[badDegreeDistri[i]] = i;
	}
	//cout << 4 << endl;
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	for (int i = 0; i < maxNodeId; i ++) 
		cout << runTimeDegree[badDegreeDistri[i]] << endl;
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
	delete [] runTimeDegree;
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
