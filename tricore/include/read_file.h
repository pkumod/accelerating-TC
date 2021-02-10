#pragma once

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "log.h"

class ReadFile {
public:
    explicit ReadFile(const char *filepath) {
        struct stat statbuf;
        if ((fdin = open(filepath, O_RDONLY)) < 0) {
            log_error("open file error");
            exit(-1);
        }
        if ((fstat(fdin, &statbuf)) < 0) {
            log_error("get file length error");
            exit(-1);
        }
        len = statbuf.st_size;
        if ((src = mmap(nullptr, len, PROT_READ, MAP_SHARED, fdin, 0)) == (void *) -1) {
            log_error("mmap file error");
            exit(-1);
        }
    }

    void release() {
        if (fdin != -1) {
            close(fdin);
            munmap(src, len);
        }
        fdin = -1;
        src = nullptr;
    }

    ~ReadFile() {
        release();
    }

    uint64_t getLen() {
        return len;
    }

    void *getAddr() {
        return src;
    }

private:
    int fdin;
    void *src;
    uint64_t len;
};
