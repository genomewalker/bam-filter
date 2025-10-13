#ifndef COVERAGE_KHASH_H
#define COVERAGE_KHASH_H

#include "khash.h"

// Define the hash map type for coverage tracking
KHASH_MAP_INIT_INT(coverage, int32_t)

// Type alias for Cython compatibility
typedef khash_t(coverage) kh_coverage_t;

#ifdef __cplusplus
extern "C" {
#endif

// Helper function for getting value pointer (needed for Cython)
static inline int* kh_val_coverage_ptr(kh_coverage_t* h, khint_t k) {
    return &kh_value(h, k);
}

#ifdef __cplusplus
}
#endif

#endif /* COVERAGE_KHASH_H */
