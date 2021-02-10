#include<iostream>
#include<string.h>
#include<fstream>
#include<stdlib.h>
#include<vector>
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
	int src, dst;
	while(fRead >> src >> dst) {
		maxNodeId = (src > maxNodeId) ? src : maxNodeId;
		maxNodeId = (dst > maxNodeId) ? dst : maxNodeId;
		allEdges.push_back(make_pair(src,dst));
	}
	fRead.close();
	
	float dAve = (float)allEdges.size()/(maxNodeId+1);

	int *lp_degree = new int [maxNodeId+1];
	int *degree = new int [maxNodeId+1];
	memset(lp_degree,0,sizeof(int)*(maxNodeId+1));
	memset(degree,0,sizeof(int)*(maxNodeId+1));
	for (int i = 0; i < allEdges.size(); i ++) {
		lp_degree[allEdges[i].first]++;
		degree[allEdges[i].first] ++;
		degree[allEdges[i].second]++;
	}
	
	float f_lp = 0;
	for (int i = 0; i <= maxNodeId; i ++) 
		f_lp += abs(dAve-lp_degree[i]);
	cout << "the lp base value is " << f_lp << endl;


	int *runTimeDegree = new int [maxNodeId+1];
	memset(runTimeDegree,0,sizeof(int)*(maxNodeId+1));
	for (int i = 0; i < allEdges.size(); i ++) {
		int src = allEdges[i].first;
		int dst = allEdges[i].second;
		if (degree[src] > degree[dst] || 
			(degree[src] == degree[dst] && src < dst)) {
	//		allEdges[i].first = dst;
	//		allEdges[i].second = src;
			runTimeDegree[dst]++;
		}
		else
			runTimeDegree[src]++;
	}
	float f_degree = 0;
	for (int i = 0; i < maxNodeId; i ++) 
		f_degree += abs(runTimeDegree[i] - dAve);
	cout << "the degree based value is " << f_degree << endl;
	return 0;
}
