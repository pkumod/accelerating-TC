#include "log.h"

static struct {
    void *udata;
    log_LockFn lock;
    FILE *fp;
    int level;
    int quiet;
    int test;
} L;

static std::mutex global_log_mutex;

static const char *level_names[] = {
        "TRACE", "DEBUG", "INFO", "WARN", "ERROR", "FATAL"
};


static void lock() {
    if (L.lock) {
        L.lock(L.udata, 1);
    }
}

static void unlock() {
    if (L.lock) {
        L.lock(L.udata, 0);
    }
}

void log_set_udata(void *udata) {
    L.udata = udata;
}

void log_set_lock(log_LockFn fn) {
    L.lock = fn;
}

void log_set_fp(FILE *fp) {
    L.fp = fp;
}

void log_set_level(int level) {
    L.level = level;
}

void log_set_quiet(int enable) {
    L.quiet = enable ? 1 : 0;
}

void log_log(int level, const char *file, int line, const char *fmt, ...) {
    if (level < L.level) {
        return;
    }
    using namespace std::chrono;
    time_point<high_resolution_clock> clock_now = high_resolution_clock::now();
    {
        std::unique_lock<std::mutex> lock_global(global_log_mutex);
        /* Acquire lock */
        lock();

        /* Get current time */
        time_t t = time(nullptr);
        struct tm *lt = localtime(&t);

        /* Log to stderr */
        if (!L.quiet) {
            va_list args;
            char buf[16];
            buf[strftime(buf, sizeof(buf), "%H:%M:%S", lt)] = '\0';
#ifdef LOG_USE_COLOR
            fprintf(
                    stderr, "%s %s%-5s\x1b[0m \x1b[90m(ts: %.6lf) %s:%d:\x1b[0m ",
                    buf, level_colors[level], level_names[level],
                    duration_cast<nanoseconds>(clock_now.time_since_epoch()).count() / 1e9, file, line);
#else
            fprintf(stderr, "%s %-5s (ts: %.6lf) %s:%d: ", buf, level_names[level],
                    duration_cast<nanoseconds>(clock_now.time_since_epoch()).count() / 1e9, file, line);
#endif
            va_start(args, fmt);
            vfprintf(stderr, fmt, args);
            va_end(args);
            fprintf(stderr, "\n");
        }

        /* Log to file */
        if (L.fp) {
            va_list args;
            char buf[32];
            buf[strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", lt)] = '\0';
            fprintf(L.fp, "%s %-5s (ts: %.6lf) %s:%d: ", buf, level_names[level],
                    duration_cast<nanoseconds>(clock_now.time_since_epoch()).count() / 1000000000.0, file, line);
            va_start(args, fmt);
            vfprintf(L.fp, fmt, args);
            va_end(args);
            fprintf(L.fp, "\n");
            fflush(L.fp);
        }

        /* Release lock */
        unlock();
    }
}