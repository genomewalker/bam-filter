#ifndef TAXONOMY_KHASH_H
#define TAXONOMY_KHASH_H

#include "khash.h"

/* Initialize khash types for taxonomy module */

/* String -> int32_t hash table for accession mapping */
KHASH_MAP_INIT_STR(str, int32_t)

/* int32_t -> int32_t hash table for LCA cache index */
KHASH_MAP_INIT_INT(int32, int32_t)

#endif /* TAXONOMY_KHASH_H */
