# cython: language_level=3

"""
Cython header file for taxonomy_db module.

This allows other .pyx files to cimport and use the taxonomy functions
directly at C speed without Python overhead.

Usage in other .pyx files:
    from bam_filter.taxonomy_db cimport (
        TaxonomyDB, AccessionMap, LCACache,
        compute_lca_nogil, query_lca_cache
    )
"""

from libc.stdint cimport int32_t, uint32_t, int64_t

# Forward declarations of C structures
cdef struct TaxNode:
    int32_t taxid
    int32_t parent_taxid
    int32_t rank_id
    int32_t name_offset
    int32_t name_length
    int32_t depth

cdef struct TaxonomyDB:
    TaxNode* nodes
    int32_t* taxid_to_idx
    char* names_buffer
    char** rank_names
    int32_t n_nodes
    int32_t n_ranks
    int32_t max_taxid
    int32_t root_idx
    int64_t names_buffer_size

cdef struct AccessionMap:
    void* acc_hash  # kh_str_t* (opaque to avoid khash dependency)
    int32_t n_entries

cdef struct LCACache:
    int32_t* lca_matrix
    int32_t* cached_taxids
    int32_t n_cached
    void* taxid_to_cache_idx  # kh_int32_t* (opaque)

# C-level functions (nogil for maximum performance)
cdef TaxonomyDB* create_taxonomy_db() nogil
cdef void free_taxonomy_db(TaxonomyDB* db) nogil
cdef AccessionMap* create_accession_map() nogil
cdef void free_accession_map(AccessionMap* amap) nogil
cdef LCACache* create_lca_cache(TaxonomyDB* db, int32_t* cached_taxids, int32_t n_cached) nogil
cdef void free_lca_cache(LCACache* cache) nogil

# Core LCA computation (nogil for parallel usage)
cdef int32_t compute_lca_nogil(TaxonomyDB* db, int32_t taxid1, int32_t taxid2) nogil
cdef int32_t query_lca_cache(LCACache* cache, int32_t taxid1, int32_t taxid2) nogil

# Lineage string generation (Greengenes-style)
cdef int build_lineage_string_nogil(TaxonomyDB* db, int32_t taxid, char* buffer, int32_t buffer_size) noexcept nogil

# Rank utilities
cdef int32_t get_rank_id(str rank_str)
cdef str get_rank_name(int32_t rank_id)

# Python wrapper classes (for cimport in other .pyx files)
cdef class TaxonomyDatabase:
    cdef TaxonomyDB* db
    cdef LCACache* lca_cache

    # Static factory method
    @staticmethod
    cdef TaxonomyDatabase _from_c_struct(TaxonomyDB* db)

    # C-level methods (accessible to other .pyx files)
    cdef int32_t _get_parent_nogil(self, int32_t taxid) nogil
    cdef int32_t _get_rank_id_nogil(self, int32_t taxid) nogil
    cdef int32_t _compute_lca_nogil(self, int32_t taxid1, int32_t taxid2) nogil

cdef class AccessionMapping:
    cdef AccessionMap* amap

    # Static factory method
    @staticmethod
    cdef AccessionMapping _from_c_struct(AccessionMap* amap)

    # C-level method for fast lookup from other .pyx files
    cdef int32_t _get_taxid_nogil(self, const char* accession) nogil
