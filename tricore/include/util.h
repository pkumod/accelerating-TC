#pragma once

#include <string>
#include <chrono>
#include <cstdio>
#include <cstdarg>
#include <vector>
#include <algorithm>

#include "log.h"

using namespace std;

#define FIRST(x) (uint32_t(x >> 32u))
#define SECOND(x) (uint32_t(x))
#define MAKEEDGE(x, y) ((uint64_t)(uint64_t(x) <<32u | uint64_t(y)))

inline string getFileName(int argc, char **argv) {
    if (argc < 3) {
        log_info("usage: %s -f graph_name", argv[0]);
        exit(-1);
    }
    return string(argv[2]);
}

template<class ForwardIt, class UnaryPredicate>
inline ForwardIt swap_if(ForwardIt first, ForwardIt last, UnaryPredicate p)
{
    if (first != last) {
        for (ForwardIt i = first; i != last; ) {
            if (p(*i)) {
                if (first != i) {
                    swap(*i, *first);
                }
                ++first;
            }
            ++i;
        }
    }
    return first;
}
