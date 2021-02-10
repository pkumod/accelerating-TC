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
struct Unit{
	int degree;
	int id;
	bool operator < (const Unit u) const{
		return (degree < u.degree)||(degree == u.degree&& id < u.id);
	}
};
using namespace std;
int main(int argc, const char *argv[]) {
	unsigned int maxNode = 0;
	ofstream fwrite(argv[2],ios::out|ios::binary);
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
	for (unsigned int i = 0; i < edgeNum; i ++) {
		Edge& e = edgeList[i];
		maxNode = e.src > maxNode ? e.src: maxNode;
		maxNode = e.dst > maxNode ? e.dst: maxNode;
	}
	unsigned nodeNum = maxNode + 1;
	cout << "the nodeNum is " << nodeNum << endl;
 	int *degree = new int [nodeNum];
	memset(degree, 0, sizeof(int)*nodeNum);
	for (unsigned int i = 0; i < edgeNum; i ++) {
		Edge& e = edgeList[i];
		degree[e.src]++;
		degree[e.dst]++;
	}
	struct Unit* units = new Unit[nodeNum];
	for (int i = 0; i < nodeNum; i ++){
		units[i].degree = degree[i];
		units[i].id = i;
	}
	cout << "start sort" <<endl;
	sort(units,units+nodeNum);
	cout << "end sort" << endl;
	/*for (int i = 0; i < 10; i ++) {
		int index = nodeNum*i/10;
		cout << "the index is " << index << endl;
		cout << units[index].degree << " " << units[index].id << endl;
	}*/
	int * lookUpTable = new int [nodeNum];
	for (int i = 0; i < nodeNum; i++) {
		lookUpTable[units[i].id] = nodeNum-i-1;
	}
	Edge e;
	for (unsigned i = 0; i < edgeNum; i ++) {
		e.src = lookUpTable[edgeList[i].src];
		e.dst = lookUpTable[edgeList[i].dst];
		fwrite.write((char *)&e,sizeof(Edge));
	}
	fwrite.close();

	delete [] lookUpTable;
	delete [] degree;
	delete [] units;
	delete [] edgeList;
	return 0;
}










