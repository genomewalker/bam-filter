// Lightweight nogil-friendly logging & timing utilities
// Header-only implementation: static inline functions so they can be used
// directly from Cython-generated C code without extra linking.
#ifndef BAM_FILTER_NOGIL_LOG_H
#define BAM_FILTER_NOGIL_LOG_H

#include <time.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

// Return monotonic time in nanoseconds
static inline unsigned long long nogil_log_now_ns(void) {
    struct timespec ts;
    // CLOCK_MONOTONIC is suitable for measuring intervals
    if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) {
        return 0ULL;
    }
    return (unsigned long long)ts.tv_sec * 1000000000ULL + (unsigned long long)ts.tv_nsec;
}

// Verbosity is read from process environment variable BAM_FILTER_VERBOSITY
// If not set, default to 1.
static inline int nogil_log_get_verbosity(void) {
    const char *s = getenv("BAM_FILTER_VERBOSITY");
    if (!s) return 1;
    return atoi(s);
}

// Set verbosity for the whole process (uses setenv so other TUs see it)
static inline void nogil_log_set_verbosity(int v) {
    char buf[16];
    int n = snprintf(buf, sizeof(buf), "%d", v);
    if (n > 0) {
        // best-effort; ignore return
        setenv("BAM_FILTER_VERBOSITY", buf, 1);
    }
}

// Simple message logger that is safe to call without the Python GIL.
// Writes a timestamp and the provided message to stderr using write().
static inline void nogil_log_msg_level(int level, const char *msg) {
    if (level > nogil_log_get_verbosity()) return;
    unsigned long long ns = nogil_log_now_ns();
    unsigned long long sec = ns / 1000000000ULL;
    unsigned long long ms = (ns / 1000000ULL) % 1000ULL;
    char buf[2048];
    int n = snprintf(buf, sizeof(buf), "[%llu.%03llu] %s\n", sec, ms, msg);
    if (n > 0) {
        // Use write to be safe in nogil sections (POSIX)
        (void)write(STDERR_FILENO, buf, (size_t)n);
    }
}

// Log elapsed time between two nanosecond timestamps with a label.
static inline void nogil_log_elapsed_level(int level, const char *label, unsigned long long start_ns, unsigned long long end_ns) {
    if (level > nogil_log_get_verbosity()) return;
    unsigned long long elapsed_us = 0ULL;
    if (end_ns >= start_ns) elapsed_us = (end_ns - start_ns) / 1000ULL;
    char buf[256];
    int n = snprintf(buf, sizeof(buf), "%s: %llu us\n", label ? label : "elapsed", elapsed_us);
    if (n > 0) {
        (void)write(STDERR_FILENO, buf, (size_t)n);
    }
}

#ifdef __cplusplus
}
#endif

#endif // BAM_FILTER_NOGIL_LOG_H
