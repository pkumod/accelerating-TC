#include<iostream>
#include<fstream>
struct Edge{
	int src;
	int dst;
};
using namespace std;
int main(int argc, const char * argv[]) {
	fstream fRead;
	fRead.open(argv[1],ios::in);
	ofstream fWrite(argv[2],ios::binary);
	char e;
	int src,dst,label;
	Edge tmpEdge;
	cout << "sizeof Edge is " << sizeof(Edge) << endl;
	while(fRead >> src >> dst){ // >> label) {
		tmpEdge.src = src;
		tmpEdge.dst = dst;
		fWrite.write((char*)(&tmpEdge),sizeof(struct Edge));
	}
	fRead.close();
	fWrite.close();
	return 0;
}
	
