#include<iostream>
#include<string.h>
#include<fstream>
using namespace std;
int main(int argc, const char * argv[]){
	fstream fRead;
	fRead.open(argv[1],ios::in);
	ofstream fWrite(argv[2],ios::out);
	int src,dst;
	int *flag = new int [60000000];
	memset(flag,-1,sizeof(int)*60000000);
	int maxValue = 0;
	int count = 0;
	while (fRead >> src >> dst) {
		if (flag[src] == -1)
			flag[src] = count++;
		if (src > maxValue) maxValue = src;
		if (dst > maxValue) maxValue = dst;
	//	cout << src << " " << dst << endl;
	}
	fRead.close();
	//for (int i = 0; i < maxValue; i++){
	//	if (flag[i] == 0)
	//		cout << "******** " << i << endl;
	//}
	cout << "the maxvalue is " << maxValue << endl;

	for (int i = 0; i <= maxValue; i++) {
		//cout << count << " " <<  flag[i] << endl;
		if (flag[i] == -1)
			flag[i] = count++;
	}
	fRead.open(argv[1],ios::in);
	while (fRead >> src >> dst) {
		fWrite << flag[src] << " " << flag[dst] << '\n';
	}
	fWrite.close();
	fRead.close();
	delete [] flag;
	return 0;
}
