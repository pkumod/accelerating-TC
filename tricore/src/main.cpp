#include <vector>
#include <cstring>

#include "log.h"
#include "util.h"
#include "read_file.h"
#include "clock.h"
#include "tricount.h"

using namespace std;

#pragma ide diagnostic ignored "openmp-use-default-none"

int main(int argc, char *argv[]) {
    if (string(argv[argc - 1]) != "debug") {
        log_set_quiet(true);
    }
    string arg_string;
    for (int i = 0; i < argc; i++) {
        arg_string += string(argv[i]) + " ";
    }
    log_info("argc: %d argv is %s", argc, arg_string.c_str());

    uint64_t gpu_mem = init_gpu();
    log_info("GPU memory: %fMB", (double) gpu_mem / 1024 / 1024);

    Clock start("Start");
    Clock read_file("ReadFile");
    Clock preprocess("Preprocess");
    Clock tri_count("TriCount");

    log_info(start.start());

    log_info(read_file.start());

    ReadFile readFile(getFileName(argc, argv).c_str());

    if (readFile.getLen() % sizeof(uint64_t)) {
        log_error("file size error: cannot be divided by 8");
        exit(-1);
    }

    uint64_t edge_count = readFile.getLen() / sizeof(uint64_t);
    log_info(read_file.count("edge_count: %llu", edge_count));

    auto edges = (uint64_t *) malloc(readFile.getLen());
    log_info(read_file.count("resize"));

    memcpy((void *) edges, readFile.getAddr(), readFile.getLen());
    readFile.release();
    log_info(read_file.count("read end"));

    log_info(preprocess.start());

    uint32_t node_count = 0;
#pragma omp parallel for reduction(max : node_count)
    for (uint64_t i = 0; i < edge_count; ++i) {
        auto &e = edges[i];
        node_count = max(node_count, FIRST(e));
        node_count = max(node_count, SECOND(e));
    }
    node_count++;
    log_info(preprocess.count("node_count: %d", node_count));

    auto deg = (uint32_t *)malloc(node_count * sizeof(uint32_t));
    cal_degree(edges, edge_count, deg, node_count);
    log_info(preprocess.count("cal degree"));

    redirect_edges(edges, edge_count, deg, node_count);
    log_info(preprocess.count("redirect_edges"));

    free(deg);
    log_info(preprocess.count("free deg"));

    uint64_t part_num = cal_part_num(edges, edge_count, node_count);
    log_info("split into part_num: %lu", part_num);

    auto *node_split = new uint32_t[part_num + 1];
    node_split[0] = 0;
    for (uint64_t i = 1; i <= part_num; i++) {
        if (i == part_num) {
            node_split[i] = node_count;
        } else {
            node_split[i] = node_count / part_num * i;
        }
    }
    auto **node_index = new uint32_t*[part_num];
    auto *adj_count = new uint64_t[part_num + 1];
    split(edges, edge_count, part_num, node_index, node_split, adj_count);
    //edges has been changed
    log_info(preprocess.count("split"));

    log_info(tri_count.start());
    auto *adj = reinterpret_cast<uint32_t *>(edges);
    auto all_sum = tricount_gpu(part_num, adj, adj_count, node_split, node_index);
    log_info(tri_count.count("triangle count: %lu", all_sum));

    log_info(start.count("all time"));
    printf("There are %lu triangles in the input graph.\n", all_sum);
    return 0;
}
