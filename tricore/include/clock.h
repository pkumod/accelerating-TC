#pragma once

#include <string>
#include <chrono>
#include <cstdio>
#include <cstdarg>

using namespace std;
using namespace std::chrono;

class Clock {
public:
    explicit Clock(const string &str) {
        s = str;
    }

    const char *start() {
        start_time = duration_cast<microseconds>(high_resolution_clock::now().time_since_epoch()).count();
        sprintf(str, "[\033[1m\033[44;30m%-12s\033[0m] Start...", s.c_str());
        return str;
    }

    const char *count(const char *fmt = "", ...) {
        va_list args;
        char str2[1000];
        va_start(args, fmt);
        vsprintf(str2, fmt, args);
        va_end(args);
        uint64_t end_time = duration_cast<microseconds>(high_resolution_clock::now().time_since_epoch()).count();
        double t = double(end_time - start_time) / 1e6;
        sprintf(str, "[\033[1m\033[44;31m%-12s\033[0m] %.6lfs   %s", s.c_str(), t, str2);
        return str;
    }

private:
    char str[1000]{};
    string s{};
    uint64_t start_time{};
};
