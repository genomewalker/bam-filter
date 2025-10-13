// Minimal nogil-safe logging and timing helpers for bam-filter
// Designed for inclusion from Cython modules in nogil sections.
// Uses clock_gettime(CLOCK_MONOTONIC) for monotonic timestamps and
// write(2) to write messages to stderr (async-signal-safe, no malloc).

#ifndef BAM_FILTER_C_LOGGING_H
#define BAM_FILTER_C_LOGGING_H

#include <time.h>
#include <unistd.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

enum {
    BF_LOG_LEVEL_QUIET = -1,
    BF_LOG_LEVEL_SUMMARY = 0,
    BF_LOG_LEVEL_INFO = 1,
    BF_LOG_LEVEL_DEBUG = 2,
    BF_LOG_LEVEL_TRACE = 3
};

// Fallback for platforms where CLOCK_MONOTONIC isn't available
#ifndef CLOCK_MONOTONIC
#define CLOCK_MONOTONIC 1
#endif

static inline double bf_monotonic_seconds(void) {
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0) {
        return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
    }
    return 0.0;
}

// Write a short log message to stderr in a nogil-safe way.
// The function composes a timestamp and the provided message into a
// fixed-size buffer on the stack and calls write(2). No allocation.
// ANSI color codes for terminal tags (disabled by default for plain logs).
// Set these to empty strings to avoid escape sequences in non-TTY logs.
#define BF_COLOR_RESET ""
#define BF_COLOR_CYAN  ""
#define BF_COLOR_YELLOW ""

static inline int bf_get_verbosity(void);
static inline int bf_should_log(int level);

static inline void bf_nogil_log(const char* tag, const char* msg) {
    char buf[512];
    const char* text = msg ? msg : "";
    size_t text_len = strlen(text);
    const char* line_end = (text_len > 0 && text[text_len - 1] == '\n') ? "" : "\n";
    int show_tag = (tag != NULL && tag[0] != '\0' && bf_get_verbosity() >= BF_LOG_LEVEL_INFO);
    int n;
    if (show_tag) {
        const char* color = BF_COLOR_CYAN;
        n = snprintf(buf, sizeof(buf), "%s%-12s%s | %s%s", color, tag, BF_COLOR_RESET, text, line_end);
    } else {
        n = snprintf(buf, sizeof(buf), "%s%s", text, line_end);
    }
    if (n > 0) {
        if (n > (int)sizeof(buf)) n = (int)sizeof(buf);
        ssize_t w = write(2, buf, (size_t)n);
        (void)w; // ignore result; safe in nogil
    }
}

// Write a formatted log to stderr. Limited formatting: takes one long and one string.
static inline void bf_nogil_log_fmt(const char* tag, const char* fmt_msg, long v) {
    char buf[512];
    int n;
    int show_tag = (tag != NULL && tag[0] != '\0' && bf_get_verbosity() >= BF_LOG_LEVEL_INFO);
    if (show_tag) {
        const char* color = BF_COLOR_YELLOW;
        n = snprintf(buf, sizeof(buf), "%s%-12s%s | ", color, tag, BF_COLOR_RESET);
    } else {
        n = 0;
    }
    if (n < 0) return;
    int m = snprintf(buf + n, sizeof(buf) - n, fmt_msg, v);
    if (m < 0) return;
    int total = n + m;
    if (total > (int)sizeof(buf)) total = (int)sizeof(buf);
    if (total > 0) {
        if (buf[total - 1] != '\n') {
            if (total < (int)sizeof(buf)) {
                buf[total++] = '\n';
            }
        }
        ssize_t w = write(2, buf, (size_t)total);
        (void)w;
    }
}

// Varargs formatted nogil logger: formats into a stack buffer then writes to stderr.
static inline void bf_nogil_logf(const char* tag, const char* fmt, ...) {
    char buf[1024];
    int n = 0;
    int show_tag = (tag != NULL && tag[0] != '\0' && bf_get_verbosity() >= BF_LOG_LEVEL_INFO);
    if (show_tag) {
        const char* color = BF_COLOR_CYAN;
        n = snprintf(buf, sizeof(buf), "%s%-12s%s | ", color, tag, BF_COLOR_RESET);
    }
    if (n < 0) return;

    va_list ap;
    va_start(ap, fmt);
    int m = vsnprintf(buf + n, sizeof(buf) - n, fmt, ap);
    va_end(ap);
    if (m < 0) return;

    int total = n + m;
    if (total > (int)sizeof(buf)) total = (int)sizeof(buf);
    if (total > 0) {
        if (buf[total - 1] != '\n') {
            if (total < (int)sizeof(buf)) {
                buf[total++] = '\n';
            }
        }
        ssize_t w = write(2, buf, (size_t)total);
        (void)w;
    }
}

// Varargs formatted nogil logger WITHOUT timestamp: formats into a stack buffer then
// writes to stderr. This is useful for modules that want to control timing output
// themselves (for example: print durations for key steps instead of per-line
// timestamps).
static inline void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) {
    char buf[1024];
    int n = 0;
    int show_tag = (tag != NULL && tag[0] != '\0' && bf_get_verbosity() >= BF_LOG_LEVEL_INFO);

    // Prepend tag without timestamp only if provided
    if (show_tag) {
        const char* color = BF_COLOR_YELLOW;
        n = snprintf(buf, sizeof(buf), "%s%-12s%s | ", color, tag, BF_COLOR_RESET);
    }
    if (n < 0) return;

    va_list ap;
    va_start(ap, fmt);
    int m = vsnprintf(buf + n, sizeof(buf) - n, fmt, ap);
    va_end(ap);
    if (m < 0) return;

    int total = n + m;
    if (total > (int)sizeof(buf)) total = (int)sizeof(buf);
    if (total > 0) {
        if (buf[total - 1] != '\n') {
            if (total < (int)sizeof(buf)) {
                buf[total++] = '\n';
            }
        }
        ssize_t w = write(2, buf, (size_t)total);
        (void)w;
    }
}

// Per-translation-unit verbosity control. This is intentionally a
// translation-unit-local static so it is safe to include from headers
// and usable from nogil contexts. It allows users to expose more
// detailed logs without allocating or locking.
static inline int *bf_verbosity_ptr(void) {
    static int v = 0;
    return &v;
}

// Set verbosity level (returns previous level)
static inline int bf_set_verbosity(int v) {
    int *p = bf_verbosity_ptr();
    int old = *p;
    *p = v;
    return old;
}

// Get current verbosity level
static inline int bf_get_verbosity(void) {
    return *bf_verbosity_ptr();
}

// Check whether the current verbosity permits emitting ``level``.
static inline int bf_should_log(int level) {
    return bf_get_verbosity() >= level;
}

// Conditionally log a pre-formatted message if the verbosity allows it.
static inline void bf_nogil_log_level(int level, const char* tag, const char* msg) {
    if (!bf_should_log(level)) {
        return;
    }
    bf_nogil_log(tag, msg);
}

// Verbose nogil logger: only emits output when current verbosity >= level.
// Uses the no-time formatter so callers can control timestamps/durations.
static inline void bf_nogil_logf_verbose(int level, const char* tag, const char* fmt, ...) {
    if (bf_get_verbosity() < level) return;
    char buf[1024];
    int n;
    int show_tag = (tag != NULL && tag[0] != '\0' && bf_get_verbosity() >= BF_LOG_LEVEL_INFO);
    if (show_tag) {
        const char* color = BF_COLOR_CYAN;
        n = snprintf(buf, sizeof(buf), "%s%-12s%s | ", color, tag, BF_COLOR_RESET);
    } else {
        n = 0;
    }
    if (n < 0) return;

    va_list ap;
    va_start(ap, fmt);
    int m = vsnprintf(buf + n, sizeof(buf) - n, fmt, ap);
    va_end(ap);
    if (m < 0) return;

    int total = n + m;
    if (total > (int)sizeof(buf)) total = (int)sizeof(buf);
    if (total > 0) {
        if (buf[total - 1] != '\n') {
            if (total < (int)sizeof(buf)) {
                buf[total++] = '\n';
            }
        }
        ssize_t w = write(2, buf, (size_t)total);
        (void)w;
    }
}

// Timing helper: compute and log a step duration (no timestamp). Start and
// end are monotonic seconds (as returned by bf_monotonic_seconds()). The
// message includes the step name and duration in seconds with millisecond
// precision.
static inline void bf_log_step_duration_notime(const char* tag, const char* step, double start, double end) {
    double dur = end - start;
    if (dur < 0.0) dur = 0.0;
    // Use a compact unit (ms) for shorter, geekier output when < 10s
    if (dur < 10.0) {
        long ms = (long)(dur * 1000.0 + 0.5);
        bf_nogil_logf_notime(tag, "%s: %ld ms\n", step ? step : "step", ms);
    } else {
        bf_nogil_logf_notime(tag, "%s took %.3f s\n", step ? step : "step", dur);
    }
}

#endif // BAM_FILTER_C_LOGGING_H
