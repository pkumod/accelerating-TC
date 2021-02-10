#include<iostream>
#include<fstream>
struct Edge{
	int src;
	int dst;
};
using namespace std;
int main(int argc, const char * argv[]) {
	fstream fRead;
	fRead.open(argv[1],ios::in|ios::binary);
	ofstream fWrite(argv[2],ios::out);
	char e;
	int src,dst,label;
	Edge tmpEdge;
	cout << "sizeof Edge is " << sizeof(Edge) << endl;
	while(fRead.read((char*)(&tmpEdge),sizeof(Edge))){//>> label) {
		fWrite << tmpEdge.src << " " << tmpEdge.dst << '\n'; 
	}
	fRead.close();
	fWrite.close();
	return 0;
}
	
