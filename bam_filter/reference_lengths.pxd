# cython: language_level=3

"""
Header file for Array-Based TSV Reference Length Module

SIMPLIFIED VERSION - NO HASH TABLES!
This allows other Cython modules to import and use the TSV reference functionality
with simple, reliable array-based lookups.

Key improvements over hash-based version:
- No hash collisions
- Simpler memory management
- More predictable performance
- Easier debugging
- No khash dependencies for the core functionality
"""

from libc.stdint cimport int32_t, int64_t

# ===============================================================================
# SIMPLIFIED DATA STRUCTURES - NO HASH TABLES
# ===============================================================================

cdef struct TSVReferenceEntry:
    char* reference_name
    int64_t reference_length

cdef struct TSVReferenceMap:
    TSVReferenceEntry* entries
    int32_t entry_count
    int32_t capacity
    bint owns_memory
    # NOTE: Removed name_to_index hash table - using simple array lookup instead

# ===============================================================================
# PRIMARY API FUNCTIONS - SIMPLE AND RELIABLE
# ===============================================================================

# Main functions for loading and using TSV reference data
cdef TSVReferenceMap* load_tsv_reference_file(const char* tsv_file_path) noexcept nogil
cdef int64_t lookup_reference_length(TSVReferenceMap* tsv_map, 
                                    const char* ref_name, 
                                    int64_t fallback_length) noexcept nogil
cdef void free_tsv_reference_map(TSVReferenceMap* tsv_map) noexcept nogil
cdef int32_t get_tsv_reference_count(TSVReferenceMap* tsv_map) noexcept nogil
cdef void print_tsv_reference_stats(TSVReferenceMap* tsv_map) noexcept nogil

# ===============================================================================
# LOW-LEVEL FUNCTIONS - FOR ADVANCED USAGE
# ===============================================================================

# Core map management
cdef TSVReferenceMap* create_tsv_reference_map() nogil
cdef void destroy_tsv_reference_map(TSVReferenceMap* tsv_map) noexcept nogil

# Entry management
cdef int add_tsv_reference_entry(TSVReferenceMap* tsv_map, 
                                const char* ref_name, 
                                int64_t ref_length) nogil

# Direct lookup (returns -1 if not found, no fallback)
cdef int64_t lookup_tsv_reference_length(TSVReferenceMap* tsv_map, 
                                        const char* ref_name) noexcept nogil

# ===============================================================================
# UTILITY FUNCTIONS
# ===============================================================================

# String and file utilities
cdef bint is_gzip_file(const char* file_path) nogil
cdef void trim_whitespace_inplace(char* str) nogil
cdef bint string_equals(const char* s1, const char* s2) nogil
