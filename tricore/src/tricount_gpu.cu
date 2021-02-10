#include "tricount.h"

#include <algorithm>
#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/copy.h>
#include <thrust/remove.h>
#include <thrust/unique.h>
#include <thrust/scan.h>
#include <device_launch_parameters.h>
#include <cuda_runtime.h>

#include "util.h"

#define CUDA_TRY(call)                                                          \
  do {                                                                          \
    cudaError_t const status = (call);                                          \
    if (cudaSuccess != status) {                                                \
      log_error("%s %s %d\n", cudaGetErrorString(status), __FILE__, __LINE__);  \
    }                                                                           \
  } while (0)

const int numBlocks = 1048576;
const int BLOCKSIZE = 512;//1024;

uint64_t gpu_mem;

uint64_t init_gpu() {
    cudaDeviceProp deviceProp{};
    CUDA_TRY(cudaGetDeviceProperties(&deviceProp, 0));
    gpu_mem = deviceProp.totalGlobalMem;
    log_info("numBlocks: %d  BLOCKSIZE: %d", numBlocks, BLOCKSIZE);
    return gpu_mem;
}

__global__ void all_degree_kernel(const uint64_t *edges, uint64_t edge_count, uint32_t *deg) {
    uint32_t blockSize = blockDim.x * gridDim.x;
    uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    for (uint64_t i = tid; i < edge_count; i += blockSize) {
        uint64_t edge = edges[i];
        auto first = FIRST(edge);
        auto second = SECOND(edge);
        atomicAdd(deg + first, 1);
        atomicAdd(deg + second, 1);
    }
}

void cal_degree(const uint64_t *edges, uint64_t edge_count, uint32_t *deg, uint32_t node_count) {
    uint64_t use_mem = node_count * sizeof(uint32_t) + 1024 * 1024 * 256;
    uint64_t edge_block = (gpu_mem - use_mem) / sizeof(uint64_t);
    uint64_t split_num = edge_count / edge_block + 1;
    edge_block = edge_count / split_num;
    uint64_t *dev_edges;
    uint32_t *dev_deg;
    CUDA_TRY(cudaMalloc((void **) &dev_edges, edge_block * sizeof(uint64_t)));
    CUDA_TRY(cudaMalloc((void **) &dev_deg, node_count * sizeof(uint32_t)));
    for (uint64_t i = 0; i < edge_count; i += edge_block) {
        uint64_t copy_size = min(edge_count - i, edge_block);
        CUDA_TRY(cudaMemcpy(dev_edges, edges + i, copy_size * sizeof(uint64_t), cudaMemcpyHostToDevice));
        all_degree_kernel<<<numBlocks, BLOCKSIZE>>> (dev_edges, copy_size, dev_deg);
        CUDA_TRY(cudaDeviceSynchronize());
    }
    CUDA_TRY(cudaMemcpy(deg, dev_deg, node_count * sizeof(uint32_t), cudaMemcpyDeviceToHost));
    CUDA_TRY(cudaFree(dev_edges));
    CUDA_TRY(cudaFree(dev_deg));
}

__global__ void redirect_edges_kernel(uint64_t *edges, uint64_t edge_count, const uint32_t *deg) {
    uint32_t blockSize = blockDim.x * gridDim.x;
    uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    for (uint64_t i = tid; i < edge_count; i += blockSize) {
        uint64_t edge = edges[i];
        auto first = FIRST(edge);
        auto second = SECOND(edge);
        if (deg[first] > deg[second] || (deg[first] == deg[second] && first > second)) {
            edges[i] = MAKEEDGE(second, first);
        }
    }
}

void redirect_edges(uint64_t *edges, uint64_t edge_count, const uint32_t *deg, uint32_t node_count) {
    uint64_t use_mem = node_count * sizeof(uint32_t) + 1024 * 1024 * 256;
    uint64_t edge_block = (gpu_mem - use_mem) / sizeof(uint64_t);
    uint64_t split_num = edge_count / edge_block + 1;
    edge_block = edge_count / split_num;
    uint64_t *dev_edges;
    uint32_t *dev_deg;
    CUDA_TRY(cudaMalloc((void **) &dev_edges, edge_block * sizeof(uint64_t)));
    CUDA_TRY(cudaMalloc((void **) &dev_deg, node_count * sizeof(uint32_t)));
    CUDA_TRY(cudaMemcpy(dev_deg, deg, node_count * sizeof(uint32_t), cudaMemcpyHostToDevice));
    for (uint64_t i = 0; i < edge_count; i += edge_block) {
        uint64_t copy_size = min(edge_count - i, edge_block);
        CUDA_TRY(cudaMemcpy(dev_edges, edges + i, copy_size * sizeof(uint64_t), cudaMemcpyHostToDevice));
        redirect_edges_kernel<<< numBlocks, BLOCKSIZE>>> (dev_edges, copy_size, dev_deg);
        CUDA_TRY(cudaDeviceSynchronize());
        CUDA_TRY(cudaMemcpy(edges + i, dev_edges, copy_size * sizeof(uint64_t), cudaMemcpyDeviceToHost));
    }
    CUDA_TRY(cudaFree(dev_edges));
    CUDA_TRY(cudaFree(dev_deg));
}

uint32_t cal_part_num(const uint64_t *edges, uint64_t edge_count, uint32_t node_count) {
    uint64_t part_num = 1;
    while (true) {
        uint64_t part_edge_count = edge_count / part_num + 1;
        uint64_t part_node_count = node_count / part_num + 1;
        uint64_t tri_use_mem = part_edge_count * sizeof(uint64_t) * 2 * 115 / 100 + (part_node_count + 1) * sizeof(uint32_t) * 2 + numBlocks * BLOCKSIZE * sizeof(uint64_t);
        if (tri_use_mem < gpu_mem) {
            break;
        }
        ++part_num;
    }
    return part_num;
}

__global__ void unzip_edges_kernel(const uint64_t *edges, uint64_t edge_count, uint32_t *edges_first, uint32_t *edges_second) {
    auto from = blockDim.x * blockIdx.x + threadIdx.x;
    auto step = gridDim.x * blockDim.x;
    for (uint64_t i = from; i < edge_count; i += step) {
        uint64_t tmp = edges[i];
        edges_first[i] = FIRST(tmp);
        edges_second[i] = SECOND(tmp);
    }
}

__global__ void node_index_construct_kernel(const uint32_t *edges_first, uint64_t edge_count, uint32_t *node_index, uint32_t node_count, uint32_t start_node) {
    auto from = blockDim.x * blockIdx.x + threadIdx.x;
    auto step = gridDim.x * blockDim.x;
    for (uint64_t i = from; i <= edge_count; i += step) {
        int64_t prev = i > 0 ? (int64_t) (edges_first[i - 1] - start_node) : -1;
        int64_t next = i < edge_count ? (int64_t) (edges_first[i] - start_node) : node_count;
        for (int64_t j = prev + 1; j <= next; ++j) {
            node_index[j] = i;
        }
    }
}

struct is_self_loop : public thrust::unary_function<uint64_t, bool> {
    __host__ __device__
    bool operator()(uint64_t x) {
        return FIRST(x) == SECOND(x);
    }
};

void split(uint64_t *edges, uint64_t edge_count, uint64_t part_num, uint32_t **node_index, const uint32_t *node_split, uint64_t *adj_count) {
    uint64_t begin = 0;
    uint64_t end;
    uint64_t adj_count_all = 0;
    auto *adj = reinterpret_cast<uint32_t *>(edges);
    uint64_t *dev_edges;
    uint32_t *dev_edges_first;
    uint32_t *dev_edges_second;
    uint32_t *dev_node_index;

    for (uint64_t i = 0; i < part_num; i++) {
        log_info("split i: %lu start", i);
        uint32_t stop_node = node_split[i + 1];
        if (i == part_num - 1) {
            end = edge_count;
        } else {
            end = swap_if(edges + begin, edges + edge_count, [&](const uint64_t edge) {
                return FIRST(edge) < stop_node;
            }) - edges;
        }
        log_info("swap_if: %d start: %lu end: %lu", i, begin, end);
        adj_count[i] = end - begin;
        uint64_t copy_size = adj_count[i] * sizeof(uint64_t);
        CUDA_TRY(cudaMalloc((void **) &dev_edges, copy_size));
        CUDA_TRY(cudaMemcpy(dev_edges, edges + begin, copy_size, cudaMemcpyHostToDevice));
        thrust::device_ptr<uint64_t> dev_ptr(dev_edges);
        thrust::sort(dev_ptr, dev_ptr + adj_count[i]);
        adj_count[i] = thrust::remove_if(dev_ptr, dev_ptr + adj_count[i], is_self_loop()) - dev_ptr;
        adj_count[i] = thrust::unique(dev_ptr, dev_ptr + adj_count[i]) - dev_ptr;

        CUDA_TRY(cudaMalloc((void **) &dev_edges_first, adj_count[i] * sizeof(uint32_t)));
        CUDA_TRY(cudaMalloc((void **) &dev_edges_second, adj_count[i] * sizeof(uint32_t)));
        unzip_edges_kernel<<<numBlocks, BLOCKSIZE>>>(dev_edges, adj_count[i], dev_edges_first, dev_edges_second);
        CUDA_TRY(cudaPeekAtLastError());
        CUDA_TRY(cudaDeviceSynchronize());
        CUDA_TRY(cudaMemcpy(adj + adj_count_all, dev_edges_second, adj_count[i] * sizeof(uint32_t), cudaMemcpyDeviceToHost));
        uint32_t node_count = node_split[i + 1] - node_split[i] + 1;
        uint32_t start_node = node_split[i];
        CUDA_TRY(cudaMalloc((void **) &dev_node_index, (node_count + 1) * sizeof(uint32_t)));
        node_index_construct_kernel<<<numBlocks, BLOCKSIZE>>>(dev_edges_first, adj_count[i], dev_node_index, node_count, start_node);
        CUDA_TRY(cudaPeekAtLastError());
        CUDA_TRY(cudaDeviceSynchronize());
        node_index[i] = (uint32_t *) malloc((node_count + 1) * sizeof(uint32_t));
        CUDA_TRY(cudaMemcpy(node_index[i], dev_node_index, (node_count + 1) * sizeof(uint32_t), cudaMemcpyDeviceToHost));
        adj_count_all += adj_count[i];

        CUDA_TRY(cudaFree(dev_edges));
        CUDA_TRY(cudaFree(dev_edges_first));
        CUDA_TRY(cudaFree(dev_edges_second));
        CUDA_TRY(cudaFree(dev_node_index));
        begin = end;
    }
    vector<uint64_t> adj_count_vec(adj_count, adj_count + part_num);
    adj_count[0] = 0;
    for (uint64_t i = 1; i <= part_num; i++) {
        adj_count[i] = adj_count[i - 1] + adj_count_vec[i - 1];
    }
}

__global__ void node_index_reconstruct_kernel(uint32_t *edges_first, const uint32_t *node_index, uint32_t node_count) {
    auto from = blockDim.x * blockIdx.x + threadIdx.x;
    auto step = gridDim.x * blockDim.x;
    for (uint64_t i = from; i < node_count; i += step) {
        for (uint64_t j = node_index[i]; j < node_index[i + 1]; ++j) {
            edges_first[j] = i;
        }
    }
}

__global__ void warp_binary_kernel(const uint32_t* __restrict__ edge_m, const uint32_t* __restrict__ node_index_m, uint32_t edge_m_count, uint32_t* __restrict__ adj_m, uint32_t start_node_n, const uint32_t* __restrict__ node_index_n, uint32_t node_index_n_count, uint32_t* __restrict__ adj_n, uint64_t *results) {
    //phase 1, partition
    uint64_t count = 0;
    __shared__ uint32_t local[BLOCKSIZE];

    uint32_t i = threadIdx.x % 32;
    uint32_t p = threadIdx.x / 32;
    for (uint32_t tid = (threadIdx.x + blockIdx.x * blockDim.x) / 32; tid < edge_m_count; tid += blockDim.x * gridDim.x / 32) {
        uint32_t node_m = edge_m[tid];
        uint32_t node_n = adj_m[tid];
        if (node_n < start_node_n || node_n >= start_node_n + node_index_n_count) {
            continue;
        }

        uint32_t degree_m = node_index_m[node_m + 1] - node_index_m[node_m];
        uint32_t degree_n = node_index_n[node_n + 1 - start_node_n] - node_index_n[node_n - start_node_n];
        uint32_t* a = adj_m + node_index_m[node_m];
        uint32_t* b = adj_n + node_index_n[node_n - start_node_n];
        if(degree_m < degree_n){
            uint32_t temp = degree_m;
            degree_m = degree_n;
            degree_n = temp;
            uint32_t *aa = a;
            a = b;
            b = aa;
        }

        //initial cache
        local[p * 32 + i] = a[i * degree_m / 32];
        __syncthreads();

        //search
        uint32_t j = i;
        while(j < degree_n){
            uint32_t X = b[j];
            uint32_t Y;
            //phase 1: cache
            int32_t bot = 0;
            int32_t top = 32;
            int32_t r;
            while(top > bot + 1){
                r = (top + bot) / 2;
                Y = local[p * 32 + r];
                if(X == Y){
                    count++;
                    bot = top + 32;
                }
                if(X < Y){
                    top = r;
                }
                if(X > Y){
                    bot = r;
                }
            }
            //phase 2
            bot = bot * degree_m / 32;
            top = top * degree_m / 32 - 1;
            while(top >= bot){
                r = (top + bot) / 2;
                Y = a[r];
                if(X == Y){
                    count++;
                }
                if(X <= Y){
                    top = r - 1;
                }
                if(X >= Y){
                    bot = r + 1;
                }
            }
            j += 32;
        }
        __syncthreads();
    }
    results[blockDim.x * blockIdx.x + threadIdx.x] = count;
}
__global__ void warp_intersection_kernel(const uint32_t* __restrict__ edge_m, 
										 const uint32_t* __restrict__ node_index_m, 
										 uint32_t edge_m_count, 
										 uint32_t* __restrict__ adj_m, 
										 uint32_t start_node_n, 
										 const uint32_t* __restrict__ node_index_n, 
										 uint32_t node_index_n_count, 
										 uint32_t* __restrict__ adj_n, 
										 uint64_t *results) {
    //phase 1, partition
    uint64_t count = 0;
    //__shared__ uint32_t local[BLOCKSIZE];

    //uint32_t i = threadIdx.x % 32;
    //uint32_t p = threadIdx.x / 32;
    for (uint32_t tid = threadIdx.x + blockIdx.x * blockDim.x; tid < edge_m_count; tid += blockDim.x * gridDim.x) {
        uint32_t node_m = edge_m[tid];
        uint32_t node_n = adj_m[tid];
        if (node_n < start_node_n || node_n >= start_node_n + node_index_n_count) {
            continue;
        }

        uint32_t degree_m = node_index_m[node_m + 1] - node_index_m[node_m];
        uint32_t degree_n = node_index_n[node_n + 1 - start_node_n] - node_index_n[node_n - start_node_n];
        uint32_t* a = adj_m + node_index_m[node_m];
        uint32_t* b = adj_n + node_index_n[node_n - start_node_n];

        //initial cache
		int i = 0, j = 0;
		while (i < degree_m && j < degree_n) {
			if (a[i] == b[j]) {
				count ++;
				i ++;
				j ++;
			} else if (a[i] < b[j]) {
				i ++;
			} else {
				j ++;
			}
		}
        //search
    }
    results[blockDim.x * blockIdx.x + threadIdx.x] = count;
}


uint64_t tricount_gpu(uint64_t part_num, uint32_t *adj, const uint64_t *adj_count, const uint32_t *node_split, uint32_t **node_index) {
    uint32_t n_result = numBlocks * BLOCKSIZE;
    uint64_t all_sum = 0;
    uint32_t *node_index_m_dev;
    uint32_t *adj_m_dev;
    uint32_t *edge_m_dev;
    uint32_t *node_index_n_dev;
    uint32_t *adj_n_dev;
    uint32_t *edge_n_dev;
    uint64_t *dev_results;
    for (uint64_t m = 0; m < part_num; m++) {
        uint32_t start_node_m = node_split[m];
        uint32_t *node_index_m = node_index[m];
        uint32_t node_index_m_count = node_split[m + 1] - node_split[m];
        if (node_index_m_count == 0) {
            continue;
        }
        uint64_t start_adj_m = adj_count[m];
        uint32_t *adj_m = adj + start_adj_m;
        uint32_t adj_count_m = adj_count[m + 1] - adj_count[m];
        CUDA_TRY(cudaMalloc((void **) &node_index_m_dev, (node_index_m_count + 1) * sizeof(uint32_t)));
        CUDA_TRY(cudaMemcpy(node_index_m_dev, node_index_m, (node_index_m_count + 1) * sizeof(uint32_t), cudaMemcpyHostToDevice));
        CUDA_TRY(cudaMalloc((void **) &adj_m_dev, adj_count_m * sizeof(uint32_t)));
        CUDA_TRY(cudaMemcpy(adj_m_dev, adj_m, adj_count_m * sizeof(uint32_t), cudaMemcpyHostToDevice));
        CUDA_TRY(cudaMalloc((void **) &edge_m_dev, adj_count_m * sizeof(uint32_t)));
        node_index_reconstruct_kernel<<<numBlocks, BLOCKSIZE>>>(edge_m_dev, node_index_m_dev, node_index_m_count);
        CUDA_TRY(cudaDeviceSynchronize());
        for (uint64_t n = m; n < part_num; n++) {
            uint32_t start_node_n = node_split[n];
            uint32_t *node_index_n = node_index[n];
            uint32_t node_index_n_count = node_split[n + 1] - node_split[n];
            uint64_t start_adj_n = adj_count[n];
            uint32_t *adj_n = adj + start_adj_n;
            uint32_t adj_count_n = adj_count[n + 1] - adj_count[n];
            CUDA_TRY(cudaMalloc((void **) &node_index_n_dev, (node_index_n_count + 1) * sizeof(uint32_t)));
            CUDA_TRY(cudaMemcpy(node_index_n_dev, node_index_n, (node_index_n_count + 1) * sizeof(uint32_t), cudaMemcpyHostToDevice));
            CUDA_TRY(cudaMalloc((void **) &adj_n_dev, adj_count_n * sizeof(uint32_t)));
            CUDA_TRY(cudaMemcpy(adj_n_dev, adj_n, adj_count_n * sizeof(uint32_t), cudaMemcpyHostToDevice));

            CUDA_TRY(cudaMalloc((void **) &dev_results, n_result * sizeof(uint64_t)));
//            log_info("tricount_gpu_edge_kernel start");
            warp_binary_kernel<<<numBlocks, BLOCKSIZE>>>(edge_m_dev, node_index_m_dev, adj_count_m, adj_m_dev, start_node_n, node_index_n_dev, node_index_n_count, adj_n_dev, dev_results);
            CUDA_TRY(cudaDeviceSynchronize());
//            log_info("tricount_gpu_edge_kernel end");
            thrust::device_ptr<uint64_t> ptr(dev_results);
            uint64_t sum = thrust::reduce(ptr, ptr + n_result);
            log_info("m: %d n: %d sum: %lu", m, n, sum);
            all_sum += sum;
            if (m != n) {
                CUDA_TRY(cudaMalloc((void **) &edge_n_dev, adj_count_n * sizeof(uint32_t)));
                node_index_reconstruct_kernel<<<numBlocks, BLOCKSIZE>>>(edge_n_dev, node_index_n_dev, node_index_n_count);
                CUDA_TRY(cudaDeviceSynchronize());
//                log_info("tricount_gpu_edge_kernel start");
                warp_binary_kernel<<<numBlocks, BLOCKSIZE>>>(edge_n_dev, node_index_n_dev, adj_count_n, adj_n_dev, start_node_m, node_index_m_dev, node_index_m_count, adj_m_dev, dev_results);
                CUDA_TRY(cudaDeviceSynchronize());
//                log_info("tricount_gpu_edge_kernel end");
                thrust::device_ptr<uint64_t> ptr_n(dev_results);
                sum = thrust::reduce(ptr_n, ptr_n + n_result);
                log_info("m: %d n: %d sum: %lu", n, m, sum);
                all_sum += sum;
                CUDA_TRY(cudaFree(edge_n_dev));
            }
            CUDA_TRY(cudaFree(node_index_n_dev));
            CUDA_TRY(cudaFree(adj_n_dev));
            CUDA_TRY(cudaFree(dev_results));
        }
        CUDA_TRY(cudaFree(node_index_m_dev));
        CUDA_TRY(cudaFree(adj_m_dev));
        CUDA_TRY(cudaFree(edge_m_dev));
    }
    return all_sum;
}
