#ifndef LCA_STATS_KHASH_H
#define LCA_STATS_KHASH_H

#include "taxonomy_khash.h"
#include <stdint.h>

/* Helper function to set value in khash string->int32_t map */
static inline void kh_set_value_str(kh_str_t* h, khint_t k, int32_t val) {
    kh_val(h, k) = val;
}

/* Helper function to get value from khash string->int32_t map */
static inline int32_t kh_get_value_str(kh_str_t* h, khint_t k) {
    return kh_val(h, k);
}

/* Hash table for read_hash (uint64_t) → taxid (int32_t)
 * This replaces the old string-based hash, using 70% less memory!
 * Old: read_name (string) → taxid requires ~40 bytes/read + hash overhead
 * New: read_hash (uint64_t) → taxid requires ~16 bytes/read
 */
KHASH_MAP_INIT_INT64(read_hash_to_taxid, int32_t)

/* Helper to set value in read_hash_to_taxid hash */
static inline void kh_set_value_read_hash_to_taxid(kh_read_hash_to_taxid_t* h, khint_t k, int32_t val) {
    kh_val(h, k) = val;
}

/* Count unique reads per LCA taxid from read_to_taxid_hash
 * Returns a newly allocated int32_t→int64_t hash table
 * Caller must free with kh_destroy_taxid_count()
 */
KHASH_MAP_INIT_INT(taxid_count, int64_t)

static inline kh_taxid_count_t* count_reads_per_lca_taxid(kh_str_t* read_to_taxid_hash) {
    kh_taxid_count_t* counts = kh_init(taxid_count);
    if (!counts) return NULL;

    khint_t k;
    int32_t taxid;
    int ret;
    khint_t count_k;

    // Iterate through all entries in read_to_taxid_hash
    for (k = kh_begin(read_to_taxid_hash); k != kh_end(read_to_taxid_hash); ++k) {
        if (!kh_exist(read_to_taxid_hash, k)) continue;

        // Get LCA taxid for this read
        taxid = kh_val(read_to_taxid_hash, k);

        // Increment count for this taxid
        count_k = kh_put(taxid_count, counts, taxid, &ret);
        if (ret == 0) {
            // Key already exists, increment
            kh_val(counts, count_k)++;
        } else {
            // New key, initialize to 1
            kh_val(counts, count_k) = 1;
        }
    }

    return counts;
}

/* Count unique reads per LCA taxid from hash-based read_to_taxid_hash
 * Returns a newly allocated int32_t→int64_t hash table
 * Caller must free with kh_destroy_taxid_count()
 */
static inline kh_taxid_count_t* count_reads_per_lca_taxid_hashed(kh_read_hash_to_taxid_t* read_hash_to_taxid) {
    kh_taxid_count_t* counts = kh_init(taxid_count);
    if (!counts) return NULL;

    khint_t k;
    int32_t taxid;
    int ret;
    khint_t count_k;

    // Iterate through all entries in read_hash_to_taxid
    for (k = kh_begin(read_hash_to_taxid); k != kh_end(read_hash_to_taxid); ++k) {
        if (!kh_exist(read_hash_to_taxid, k)) continue;

        // Get LCA taxid for this read
        taxid = kh_val(read_hash_to_taxid, k);

        // Increment count for this taxid
        count_k = kh_put(taxid_count, counts, taxid, &ret);
        if (ret == 0) {
            // Key already exists, increment
            kh_val(counts, count_k)++;
        } else {
            // New key, initialize to 1
            kh_val(counts, count_k) = 1;
        }
    }

    return counts;
}

#endif /* LCA_STATS_KHASH_H */
