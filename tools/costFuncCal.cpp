#include<iostream>
#include<string.h>
#include<fstream>
#include<vector>
#include<iterator>
#include<algorithm>
struct Edge{
	int src;
	int dst;
};
float myabs(float a) {
	return a > 0? a:-a;
}
using namespace std;
int main(int argc, const char *argv[]) {
	ifstream fin(argv[1],ios::in|ios::binary);
	fin.seekg(0,ios::end);
	streampos Size = fin.tellg();
	fin.close();
	long int size = Size;
	unsigned edgeNum = size/sizeof(Edge);
	Edge * edgeList = new Edge[edgeNum];
	ifstream fread(argv[1],std::ios::binary);
	fread.read((char*)edgeList,sizeof(Edge)*edgeNum);
	fread.close();
	unsigned maxNode = 0;
	for (unsigned int i = 0; i < edgeNum; i ++) {
		Edge& e = edgeList[i];
		maxNode = e.src > maxNode ? e.src: maxNode;
		maxNode = e.dst > maxNode ? e.dst: maxNode;
	}
	unsigned nodeNum = maxNode + 1;
	cout << "the nodeNum is " << nodeNum << endl;
	float avgDegree = (float)edgeNum/nodeNum;
	cout << "the avgDegree is " << avgDegree << endl;
 	int *degree = new int [nodeNum];
	int *idBasedDegree = new int[nodeNum];
	int *degreeBasedDegree = new int [nodeNum];
	int *lpBasedDegree = new int [nodeNum];
	memset(degree, 0, sizeof(int)*nodeNum);
	memset(idBasedDegree, 0, sizeof(int)*nodeNum);
	memset(degreeBasedDegree,0,sizeof(int)*nodeNum);
	memset(lpBasedDegree,0,sizeof(int)*nodeNum);
	for (unsigned int i = 0; i < edgeNum; i ++) {
		Edge& e = edgeList[i];
		degree[e.src]++;
		degree[e.dst]++;
		lpBasedDegree[e.src]++;
	}
	for (unsigned int i = 0; i < edgeNum; i ++) {
		int src = edgeList[i].src;
		int dst = edgeList[i].dst;
		if (src < dst)
			idBasedDegree[src] ++;
		else
			idBasedDegree[dst] ++;
		if (degree[src] < degree[dst] ||(degree[src] == degree[dst] && src < dst))
			degreeBasedDegree[src]++;
		else
			degreeBasedDegree[dst]++;
	}
	
	float idCost = 0, idPartCost = 0;
	float degreeCost = 0, degreePartCost = 0;
	float lpCost = 0, lpPartCost = 0;
	float countThreshold = atoi(argv[2])*avgDegree;
	for (int i = 0; i < nodeNum; i ++) {
		float curValue = myabs(idBasedDegree[i]-avgDegree);
		idCost += curValue;
		if (idBasedDegree[i]-avgDegree > countThreshold)
			idPartCost += curValue;

		curValue = myabs(degreeBasedDegree[i]-avgDegree);
		degreeCost += curValue;
		if (degreeBasedDegree[i] - avgDegree > countThreshold)
			degreePartCost += curValue;

		curValue = myabs(lpBasedDegree[i]-avgDegree);
		lpCost += curValue;
		if (lpBasedDegree[i]-avgDegree > countThreshold)
			lpPartCost += curValue;
	}
	//cout << "idCost: "<< idCost << endl;
	cout << "idPartCost: " << idPartCost << endl;
	//cout << "degreeCost: " << degreeCost << endl;
	cout << "degreePartCost: " << degreePartCost << endl;
	//cout << "lpCost: " << lpCost << endl;
	cout << "lpPartCost: " << lpPartCost << endl;
	delete [] degree;
	delete [] degreeBasedDegree;
	delete [] idBasedDegree;
	delete [] lpBasedDegree;
	delete [] edgeList;
	return 0;
}










