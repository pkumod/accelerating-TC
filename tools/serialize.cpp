#include<iostream>
#include<string.h>
#include<fstream>
using namespace std;
struct edge {
	unsigned int s;
	unsigned int d;
};
int main(int argc, const char * argv[]){
	ifstream fRead(argv[1],ios::binary);
	ofstream fWrite(argv[2],ios::binary);
	unsigned int src,dst;
	unsigned int *flag = new unsigned int [200000000];
	memset(flag,0,sizeof(unsigned int)*200000000);
	unsigned int maxValue = 0;
	struct edge E;
	while (fRead.read((char *)&E, sizeof(edge))) {
		src = E.s;
		dst = E.d;
		flag[src] = 1;
		flag[dst] = 1;
		if (src > maxValue) maxValue = src;
		if (dst > maxValue) maxValue = dst;
		//cout << src << " " << dst << endl;
	}
	fRead.close();
	/*for (int i = 0; i < maxValue; i++){
		if (flag[i] == 0)
			cout << "******** " << i << endl;
	}*/
	cout << "the maxvalue is " << maxValue << endl;
	unsigned int count = 0;
	for (unsigned int i = 0; i <= maxValue; i++) {
		//cout << count << " " <<  flag[i] << endl;
		unsigned int curFlag = flag[i];
		flag[i] = count;
		if (curFlag == 0)
			count ++;
	}
	fRead.open(argv[1],ios::binary);
	while (fRead.read((char *)&E, sizeof(edge))) {
		E.s = E.s - flag[E.s];
		E.d = E.d - flag[E.d];
		fWrite.write((char *)&E,sizeof(edge));
//		fWrite << src-flag[src] << " " << dst-flag[dst] << '\n';
	}
	fWrite.close();
	fRead.close();
	delete [] flag;
	return 0;
}
