#include <iostream>
#include <fstream>
#include "cuda_runtime.h"
#include <thrust/sort.h>
#include <thrust/sequence.h>
#include <thrust/execution_policy.h>
#include <time.h>
using namespace std;
typedef unsigned long edge;
#define FIRST(x) (uint32_t(x>>32u))
#define SECOND(x) (uint32_t(x))
int main(int argc, const char *argv []) {
	ifstream fread(argv[1], std::ios::binary);
	fread.seekg(0,ios::end);
	streampos Size = fread.tellg();
	fread.close();
	uint32_t edge_num = Size/sizeof(edge);
	edge *Edge = new edge [edge_num];

	fread.open(argv[1],std::ios::binary);
	fread.read((char *)Edge, sizeof(edge)*edge_num);
	fread.close();
	cout << "the edge num is " << edge_num << endl;

	uint32_t max_node = 0;
	for (uint32_t i = 0; i < edge_num; i ++) {
		max_node = max(max_node, FIRST(Edge[i]));
		max_node = max(max_node, SECOND(Edge[i]));
	}
	uint32_t node_num = max_node + 1;
	cout << "the node num is " << node_num << endl;

	uint32_t *degree = new uint32_t [node_num];
	memset(degree, 0, sizeof(uint32_t)*node_num);
	for (uint32_t i = 0; i < node_num; i ++) {
		degree[FIRST(Edge[i])] ++;
		degree[SECOND(Edge[i])] ++;
	}
	uint32_t max_degree = 0;
	for (uint32_t i = 0; i < node_num; i ++) {
		if (degree[i] > max_degree)
			max_degree = degree[i];
	}
	cout << "the maxdegree is " << max_degree << endl;
	//build csr
	uint32_t *offset = new uint32_t [node_num + 1];
	thrust::exclusive_scan(thrust::host, degree, degree + node_num, offset);
	offset[node_num] = degree[node_num-1] + offset[node_num-1];
	uint32_t *values = new uint32_t [edge_num * 2];
	uint32_t *runtime_offset = new uint32_t [node_num];
	memcpy(runtime_offset, offset, sizeof(uint32_t)*node_num);
	for (uint32_t i = 0; i < edge_num; i ++) {
		uint32_t src = FIRST(Edge[i]);
		uint32_t dst = SECOND(Edge[i]);
		values[runtime_offset[src]++] = dst;
		values[runtime_offset[dst]++] = src;
	}
	delete [] runtime_offset;
	cout << "csr build " << endl;
	clock_t _start = clock();

	cudaEvent_t d_start,d_stop;
	cudaEventCreate(&d_start);cudaEventCreate(&d_stop);
	cudaEventRecord(d_start, 0);

	uint32_t *d_degree;
	uint32_t *d_D;
	cudaMalloc((void **)&d_degree, sizeof(uint32_t)*node_num);
	cudaMalloc((void **)&d_D, sizeof(uint32_t)*node_num);

	thrust::sequence(thrust::device, d_D, d_D+node_num);
	cudaMemcpy(d_degree,degree,sizeof(uint32_t)*node_num,cudaMemcpyHostToDevice);
	//prepare basic vars for the implementation
	uint32_t *b = new uint32_t [max_degree + 2];
	uint32_t *D = new uint32_t [node_num];
	uint32_t *p = new uint32_t [node_num];
	uint32_t *degree_tmp = new uint32_t [node_num];
	//memcpy(degree_tmp, degree, sizeof(uint32_t)*node_num);
	thrust::sort_by_key(thrust::device, d_degree, d_degree + node_num, d_D);
	cudaMemcpy(D, d_D, sizeof(uint32_t)*node_num, cudaMemcpyDeviceToHost);
	clock_t to_sort = clock();
	cout << "end sort uses " << (double)(to_sort - _start)/CLOCKS_PER_SEC << " s." << endl;

	cudaEventRecord(d_stop,0);cudaEventSynchronize(d_stop);
	float elapsedTime = 0;
	cudaEventElapsedTime(&elapsedTime,d_start,d_stop);
	cout << "sort uses " << elapsedTime << " ms in GPU time clock" << endl;
	cudaFree(d_D);
	cudaFree(d_degree);

	b[0] = 0;
	uint32_t last_degree = 0;
	for (uint32_t i = 0; i < node_num; i ++) {
		uint32_t local_node = D[i];
		p[local_node] = i;
		if (degree[local_node] > last_degree) {
			for (uint32_t j = last_degree + 1; j <= degree[local_node]; j ++)
				b[j] = i;
			last_degree = degree[local_node];
		}
	}
	b[max_degree + 1] = node_num;
	memcpy(degree_tmp, degree, sizeof(uint32_t)*node_num);
	clock_t prepared = clock();
	cout << "vars prepared, and uses " << (double)(prepared - _start)/CLOCKS_PER_SEC << " s." << endl;
		clock_t start = clock();
	//fa algorithm
	for (uint32_t i = 0; i < node_num; i ++) {
		uint32_t v = D[i];
		for (uint32_t j = 0; j < degree_tmp[v];  j++) {
			uint32_t u = values[offset[v] + j];
			if (degree[u] > degree[v]) {
				uint32_t du = degree[u];
				uint32_t pu = p[u];
				uint32_t pw = b[du];
				uint32_t w = D[pw];
				if (u != w) {
					D[pu] = w; 
					D[pw] = u; 
					p[u] = pw; 
					p[w] = pu;
				}
				b[du]++; degree[u]--;
			}
		}
	}
	cout << "the algorithm using flatten array is finished" << endl;
	clock_t end = clock();
	cout << "the algorithm uses " << (double)(end-start)/CLOCKS_PER_SEC << " s." << endl;
	cout << "the total process uses " << (double)(end-_start)/CLOCKS_PER_SEC << " s." << endl;

	edge e;
	ofstream fwrite(argv[2], std::ios::binary);
	for (uint32_t i = 0; i < edge_num; i ++) {
		uint32_t src = FIRST(Edge[i]);
		uint32_t dst = SECOND(Edge[i]);
		if (degree[src] > degree[dst] ||
			(degree[src] == degree[dst] && src > dst)) {
			uint32_t tmp = src;
			src = dst;
			dst = tmp;
		}
		e = ((edge)src << 32) | dst;
		fwrite.write((char *)&e, sizeof(edge));
	}
	fwrite.close();

	delete [] values;
	delete [] offset;
	delete [] Edge;
	delete [] b;
	delete [] D;
	delete [] p;
	delete [] degree;
	delete [] degree_tmp;
}
