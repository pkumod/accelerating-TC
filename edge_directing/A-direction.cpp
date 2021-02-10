#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/plod_generator.hpp>
#include <boost/random/linear_congruential.hpp>

#include <iostream> // for std::cout
#include <utility> // for std::pair
#include <algorithm> // for std::for_each
#include <boost/graph/graph_traits.hpp>

#include<fstream>
#include<set>
#include<algorithm>
#include<vector>
#include<cstdlib>
#include<queue>
#include<map>
#include<iterator>
#include<ctime>
#include<omp.h>
//#include<pair>
#define EDGELIST
#define EDGEBINARY
//#define DEGREE_EDGEBINARY
//#define DEGREE_EDGELIST

#define READFROMFILE
//#define MAKEPOWERLAWGRAPH
using namespace boost;
//using namespace std;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
typedef boost::plod_iterator<boost::minstd_rand, Graph> SFGen;

float optTime;
int thread_num;
double thre_param;
int * peel(std::vector<std::pair<int,int> > & edges, unsigned int edgeNum, int nodeNum) {
	
	float threshold = 2*((float)edgeNum)/nodeNum;
	int * priority = new int [nodeNum];
	memset(priority, -1, sizeof(int)*nodeNum);
	int * degree = new int [nodeNum];
	memset(degree, 0, sizeof(int)*nodeNum);
	for (unsigned i = 0; i < edgeNum; i ++) {
		degree[edges[i].first]++;
		degree[edges[i].second]++;
	}
	int *runTimeDegree = new int [nodeNum];
	memcpy(runTimeDegree,degree,sizeof(int)*nodeNum);
	long long  * edgeOffsetPtr = new long long [nodeNum];
	edgeOffsetPtr[0] = 0;
	for (int i = 0; i < nodeNum-1; i++) {
		edgeOffsetPtr[i+1] = edgeOffsetPtr[i]+degree[i];
	}
	std::cout << "*******" << std::endl;
	std::cout << edgeNum << std::endl;
	//std::cout << (((long long)edgeNum)<<1) << std::endl;
	//while(1){}
	int * allEdges = new int[((long long)edgeNum)*2];
	long long  * edgeCurPos = new long long  [nodeNum];
	memcpy(edgeCurPos, edgeOffsetPtr, sizeof(long long)*nodeNum);
//	for ( int i = 0; i < nodeNum; i ++) {
//		edgeCurPos[i] = edgeOffsetPtr[i];
//	}
	std::cout << "vars prepared" << std::endl;
	for (unsigned int i = 0; i < edgeNum; i++) {
		int src = edges[i].first;
		int dst = edges[i].second;
		allEdges[edgeCurPos[src]] = dst;
		edgeCurPos[src] ++;
		allEdges[edgeCurPos[dst]] = src;
		edgeCurPos[dst] ++;
	}
	std::cout << "all edges duplicated" << std::endl;
	delete [] edgeCurPos;
	int priorityNum = 0;
	int threshold_old = threshold;
	int peeledOffNodes = 0;
	int *Frontier = new int[nodeNum];
	int *nextFrontier = new int [nodeNum];
	memset(Frontier, 0,sizeof(int)*nodeNum);
	memset(nextFrontier,0,sizeof(int)*nodeNum);
	int fIndex = 0,nIndex = 0;
	omp_set_nested(true);
	for (int t = 0; ; t ++) {
		//clock_t startTime = clock();
#pragma omp parallel for num_threads(2*thread_num-1)
		for (int i = 0; i < nodeNum; i ++) {
			if (runTimeDegree[i] <= threshold && priority[i] == -1) {
				int curIndex = __sync_add_and_fetch(&fIndex,1);
				Frontier[curIndex] = i;
//#pragma omp atomic
//				peeledOffNodes ++;
				priority[i] = priorityNum;
			}
		}
		double startTime = omp_get_wtime();
		priorityNum++;
		int innerLoop = 0;
		while (fIndex > nodeNum/5000 && fIndex > 200 && innerLoop <= 5) {
#pragma omp parallel for num_threads(2*thread_num-1) 
			for (int i = 0; i < fIndex; i++) {
				int curNode = Frontier[i];
				int * curOffset = allEdges + edgeOffsetPtr[curNode];
//#pragma omp parallel for num_threads(2*thread_num-1) 
				for (int j = 0; j < degree[curNode]; j ++) {
					int node = curOffset[j];
					if (priority[node] != -1) {
						if (priority[node] == priority[curNode] && 
							(runTimeDegree[node] > runTimeDegree[curNode] || (runTimeDegree[node] == runTimeDegree[curNode]&& node > curNode))) {
					   		runTimeDegree[node]--;
						}
					}else {
						runTimeDegree[node] --;
						if (runTimeDegree[node] <= threshold) {
							int curIndex = __sync_add_and_fetch(&nIndex,1);
							nextFrontier[curIndex] = node;
//#pragma omp atomic
//							peeledOffNodes++;
							priority[node] = priorityNum;
						}
					}
				}
			}
			int * tmpPtr = Frontier;
			Frontier = nextFrontier;
			nextFrontier = tmpPtr;
			fIndex = nIndex;
			nIndex = 0;
			priorityNum++;
			memset(nextFrontier,0,sizeof(int)*nodeNum);
			std::cout << "##### one peel done!####### " << fIndex << std::endl;
			innerLoop ++;
		}
		double endTime = omp_get_wtime();
		optTime += (endTime-startTime);
		std::cout << "the peel time is " << optTime << " s." << std::endl;
		bool peelDone = true;
#pragma omp parallel for 
		for (int i = 0; i < nodeNum; i ++)
			if (priority[i] == -1)
				peelDone = false;
		if (peelDone)
		//if (peeledOffNodes >= nodeNum)
			break;
		threshold += (double)threshold_old*thre_param/5;
	}
	delete [] nextFrontier;
	delete [] Frontier;

	int peelOffNodeNum = 0;
	int peeledEdges = 0;

	std::cout << "peel done, and the peeled node num is " << peelOffNodeNum << std::endl;
	
	delete [] allEdges;
	delete [] edgeOffsetPtr;
	delete [] degree;
	//delete [] priority;
	delete [] runTimeDegree;
	return priority;
}
struct node{
	int id;
	int maxOutDegreePossible;
};
struct cmp{
	bool operator()(std::pair<int,int> a, std::pair<int,int> b) {
		return a.second > b.second;
	}
};
bool cmpPair(std::pair<int,int> a, std::pair<int,int> b) {
	return b.second > a.second;
}
struct Edge{
	int src;
	int dst;
};
int main(int argc, const char *argv[])
{
	optTime = 0;
	thread_num = omp_get_num_procs();
	std::cout << "the max thread num is " << thread_num << std::endl;
	//******************************************graph generation
	int nodeNum = 0;
#ifdef MAKEPOWERLAWGRAPH
	boost::minstd_rand gen;
	// Create graph with 100 nodes 
	nodeNum = 5000000;
	Graph g(SFGen(gen, nodeNum, 1.1, 50000000), SFGen(), nodeNum);//1.4=329 1.6=195
															  //1.8=114 1.7=169

	boost::property_map<Graph, vertex_index_t>::type vertex_id = get(vertex_index, g);
#endif
#ifdef READFROMFILE
	std::fstream fRead;
	Edge curEdge;
	fRead.open(argv[1],std::ios::in|std::ios::binary);
	thre_param = (double)atoi(argv[2]);
	while(fRead.read((char*)&curEdge,sizeof(Edge))) {
		if (curEdge.src > nodeNum) nodeNum = curEdge.src;
		if (curEdge.dst > nodeNum) nodeNum = curEdge.dst;
	}
	nodeNum ++;
	fRead.close();
	std::cout << "the max nodeNum is " << nodeNum << std::endl;
	fRead.open(argv[1],std::ios::in|std::ios::binary);
#endif 
	int totalEdge = 0;

	std::vector<int> * neighbors = new std::vector<int> [nodeNum];
	//*******************************************read in graph
#ifdef MAKEPOWERLAWGRAPH
	graph_traits<Graph>::edge_iterator ei, ei_end;
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
		int s = get(vertex_id, source(*ei, g)),
 			e = get(vertex_id, target(*ei, g)) ;
		if (s == e)
			continue;
//		std::cout << "(" <<	s << "," << e << ") " << std::endl;
//		totalEdge ++;
		neighbors[s].push_back(e);
		neighbors[e].push_back(s);
	}
#endif
#ifdef READFROMFILE 
	while(fRead.read((char*)&curEdge,sizeof(Edge))) {
		int s = curEdge.src, e = curEdge.dst;
		if (s == e)
			continue;
		neighbors[s].push_back(e);
		neighbors[e].push_back(s);
	}
	fRead.close();
#endif
#pragma omp parallel for
	for (int i = 0; i < nodeNum; i ++) {
		sort(neighbors[i].begin(),neighbors[i].end());
		neighbors[i].erase(unique(neighbors[i].begin(),neighbors[i].end()),neighbors[i].end());
#pragma omp atomic
		totalEdge += neighbors[i].size();
	}
	std::cout << "the total edges is " << totalEdge/2 << std::endl;
	//for (int i = 0; i < nodeNum; i++) {
	//	for (int j = 0; j < neighbors[i].size(); j ++)
	//		std::cout << i << " " << neighbors[i][j] << std::endl;
	//}
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<peel off graph
#ifdef EDGELIST
	std::ofstream outputLPGraph("lpBasedGraph",std::ios::out);
#endif
#ifdef EDGEBINARY
	std::ofstream outputLPGraphBIN("lpBasedGraph.bin",std::ios::out|std::ios::binary);
#endif 
	std::vector<std::pair<int,int> > originEdges;
	originEdges.clear();
	for(int i = 0; i < nodeNum; i ++) {
		for (int j = 0; j < neighbors[i].size(); j ++) {
			originEdges.push_back(std::make_pair(i,neighbors[i][j]));
		}
	}
	float originThreshold = (float)totalEdge/nodeNum;
	long outputEdges = 0;
	//peel
	int *priorities = peel(originEdges,totalEdge,nodeNum);
	//output graph
	for (int i = 0; i < nodeNum; i ++) {
	//	if (priorities[i] == -1)
	//		continue;
		for (int j = 0; j < neighbors[i].size(); j++) {
			int curNgb = neighbors[i][j];
			if (//priorities[curNgb] == -1 || 
				priorities[curNgb] > priorities[i] ||
				(priorities[curNgb] == priorities[i] && i < curNgb)) {
#ifdef EDGELIST
				outputLPGraph << i << " " << curNgb << '\n';
#endif
#ifdef EDGEBINARY
				outputLPGraphBIN.write((char*)&i,sizeof(int));
				outputLPGraphBIN.write((char*)&curNgb,sizeof(int));
#endif
				outputEdges++;
			}

		}
	}
	originEdges.clear();
	for (int i = 0; i < nodeNum; i ++)
		neighbors[i].clear();
	delete [] neighbors;
	delete [] priorities;
	std::cout << "the total optTime is " << optTime << " s." << std::endl;
	std::cout << "we have dealed with " << outputEdges << " edges up to now" << std::endl;
#ifdef EDGELIST
	outputLPGraph.close();
#endif
#ifdef EDGEBINARY
	outputLPGraphBIN.close();
#endif
    return 0;
}














