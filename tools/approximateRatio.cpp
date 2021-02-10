#include<iostream>
#include<string.h>
#include<fstream>
#include<stdlib.h>
#include<vector>
#include<algorithm>
//#include<random>
#include<time.h>
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
	
	double avg_outd = (double)allEdges.size()/(maxNodeId+1);
	int dmax = 0;
	for (int i = 0; i <= maxNodeId; i ++)
			if (degree[i] > dmax)
					dmax = degree[i];
	int *percent = new int [sizeof(int)*(dmax+1)];
	memset(percent,0,sizeof(int)*(dmax+1));
	double v_low_node_size = 0, v_high_node_size = 0;
	double sumDegreeVLow = 0, sumDegreeVHigh = 0;
	for (int i = 0; i <= maxNodeId; i ++) {
			if (degree[i] <= avg_outd) {
					v_low_node_size ++;
					sumDegreeVLow += degree[i];
			} else {
					v_high_node_size ++;
					sumDegreeVHigh += degree[i];
			}
			percent[degree[i]] ++;
	}
	double sumDegree = sumDegreeVLow + sumDegreeVHigh;

	double lb_opt = 0;
	for(int d = 1; d <= avg_outd; d++)
			lb_opt += percent[d] * (avg_outd - d);		
	double lb_opt2 = 0;
	double avg_outd_vhigh = sumDegreeVHigh / v_high_node_size / 2;
	if(avg_outd_vhigh >= avg_outd)
	{
			if((sumDegreeVHigh - sumDegreeVLow) / 2 - avg_outd * v_high_node_size >= 0)
					lb_opt2 = (sumDegreeVHigh - sumDegreeVLow) / 2 - avg_outd * v_high_node_size; // Case 2(a)
			else
					lb_opt2 = 0; // Case 2(b)
	}
	else
			lb_opt2 = avg_outd * v_high_node_size - sumDegreeVHigh / 2; // Case 1

	double ub_alg_minus_opt = 0;		
	double remainEdge = 0;
	for(int d = (int)avg_outd + 1; d <= dmax; d++)
			remainEdge += percent[d] * d / 2;
	int dx = (int)avg_outd + 1;
	while(remainEdge > 0)
	{
			if(remainEdge > dx * percent[dx])
			{
					ub_alg_minus_opt += avg_outd * percent[dx];
					remainEdge -= dx * percent[dx];
					dx ++;
			}
			else
			{
					ub_alg_minus_opt += avg_outd * remainEdge / dx;
					break;
			}
	}

	double approx_ratio = 1 + ub_alg_minus_opt / (lb_opt + lb_opt2);
	cout << "The approximation ratio is " << approx_ratio << endl;
	delete [] degree;
	delete [] percent;
	allEdges.clear();
	return 0;
}
