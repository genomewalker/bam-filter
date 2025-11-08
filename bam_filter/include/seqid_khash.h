#ifndef SEQID_KHASH_H
#define SEQID_KHASH_H
#include "khash.h"
#include <stdlib.h>
#include <string.h>

/* Map types */
KHASH_MAP_INIT_INT64(seqid_map, int)
KHASH_MAP_INIT_INT64(seqid_name_map, char *)
KHASH_MAP_INIT_STR(str_map, int)
KHASH_MAP_INIT_INT64(pos_to_idx, int)
/* Partner mask map: global map type for partner->uint64_t* masks */
KHASH_MAP_INIT_INT(partner_mask, uint64_t *)
/* Edge map for leiden graph: edge_key (uint64_t) -> weight (uint32_t) */
KHASH_MAP_INIT_INT64(edge_map, unsigned int)

/* Existing size/exists macros */
#define kh_val_pos_to_idx(h, k) kh_val(h, k)
#define kh_end_seqid_map(h) kh_end(h)
#define kh_exist_seqid_map(h, k) kh_exist(h, k)
#define kh_size_seqid_map(h) kh_size(h)
#define kh_end_seqid_name_map(h) kh_end(h)
#define kh_exist_seqid_name_map(h, k) kh_exist(h, k)
#define kh_size_seqid_name_map(h) kh_size(h)
#define kh_end_str_map(h) kh_end(h)
#define kh_size_str_map(h) kh_size(h)

/* pos_to_idx macros */
#define kh_size_pos_to_idx(h) kh_size(h)
#define kh_end_pos_to_idx(h) kh_end(h)
#define kh_exist_pos_to_idx(h, k) kh_exist(h, k)

/* ADD THESE MISSING MACROS FOR pos_to_idx */
#define kh_resize_pos_to_idx(h, s) kh_resize(pos_to_idx, h, s)
#define kh_init_pos_to_idx() kh_init(pos_to_idx)
#define kh_destroy_pos_to_idx(h) kh_destroy(pos_to_idx, h)
#define kh_clear_pos_to_idx(h) kh_clear(pos_to_idx, h)
#define kh_get_pos_to_idx(h, k) kh_get(pos_to_idx, h, k)
#define kh_put_pos_to_idx(h, k, r) kh_put(pos_to_idx, h, k, r)
#define kh_del_pos_to_idx(h, k) kh_del(pos_to_idx, h, k)

#ifdef __cplusplus
extern "C" {
#endif

/* Existing value wrappers */
static inline int *kh_val_seqid_map_wrap(kh_seqid_map_t *h, khint_t k) {
  return &kh_val(h, k);
}
static inline char **kh_val_seqid_name_map_wrap(kh_seqid_name_map_t *h,
                                                khint_t k) {
  return &kh_val(h, k);
}
static inline int *kh_val_str_map_wrap(kh_str_map_t *h, khint_t k) {
  return &kh_val(h, k);
}

/* --- NEW key access + delete wrappers --- */
static inline char *kh_key_str_map_wrap(kh_str_map_t *h, khint_t k) {
  return (char *)kh_key(h, k);
}
static inline void kh_set_key_str_map_wrap(kh_str_map_t *h, khint_t k,
                                           char *new_key) {
  kh_key(h, k) = new_key;
}
static inline void kh_del_str_map_wrap(kh_str_map_t *h, khint_t k) {
  kh_del_str_map(h, k);
}

static inline int *kh_val_pos_to_idx_wrap(kh_pos_to_idx_t *h, khint_t k) {
  return &kh_val_pos_to_idx(h, k);
}
static inline void kh_destroy_pos_to_idx_wrap(kh_pos_to_idx_t *h) {
  kh_destroy_pos_to_idx(h);
}

/* Partner-mask wrappers: provide small helpers so Cython-generated code can
   call kh_val_partner_mask_wrap, kh_end_partner_mask, and
   kh_exist_partner_mask. */
static inline uint64_t **kh_val_partner_mask_wrap(kh_partner_mask_t *h,
                                                  khint_t k) {
  return &kh_val(h, k);
}

/* Map partner_mask helpers to generic khash macros */
#define kh_end_partner_mask(h) kh_end(h)
#define kh_exist_partner_mask(h, k) kh_exist(h, k)

/* Edge map macros for leiden graph */
#define kh_init_edge_map() kh_init(edge_map)
#define kh_destroy_edge_map(h) kh_destroy(edge_map, h)
#define kh_get_edge_map(h, k) kh_get(edge_map, h, k)
#define kh_put_edge_map(h, k, r) kh_put(edge_map, h, k, r)
#define kh_end_edge_map(h) kh_end(h)
#define kh_exist_edge_map(h, k) kh_exist(h, k)
#define kh_size_edge_map(h) kh_size(h)
#define kh_key_edge_map(h, k) kh_key(h, k)
#define kh_value_edge_map(h, k) kh_value(h, k)

/* Edge map value wrapper for Cython */
static inline unsigned int *kh_val_edge_map_wrap(kh_edge_map_t *h, khint_t k) {
  return &kh_val(h, k);
}

#ifdef __cplusplus
}
#endif
#endif
