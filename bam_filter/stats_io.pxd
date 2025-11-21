# cython: language_level=3
"""
Statistics output writer - C API declarations.
"""

from libc.stdio cimport FILE
from bam_filter.stats cimport RefStats, FilterConditions
from bam_filter.generic_filters cimport GenericFilters
from bam_filter.processor_types cimport sam_hdr_t

cdef extern from "zlib.h":
    ctypedef void* gzFile

# Main API functions
cdef int write_output_files_complete(
    const char* output_c,
    const char* filtered_output_c,
    RefStats* global_ref_stats,
    sam_hdr_t* header,
    int n_refs,
    FilterConditions* filters,
    GenericFilters* gfilters
) except -1 nogil

cdef int write_output_files_complete_fast(
    const char* output_c,
    const char* filtered_output_c,
    RefStats* global_ref_stats,
    sam_hdr_t* header,
    int n_refs,
    FilterConditions* filters,
    GenericFilters* gfilters,
    int compression_level,
    int compression_threads
) except -1 nogil

# Legacy compatibility functions (deprecated)
cdef int write_stats_to_file(
    FILE* fp,
    RefStats* global_ref_stats,
    sam_hdr_t* header,
    int n_refs,
    FilterConditions* filters,
    GenericFilters* gfilters,
    bint apply_filters
) except -1 nogil

cdef int write_stats_to_gzfile(
    gzFile gzfp,
    RefStats* global_ref_stats,
    sam_hdr_t* header,
    int n_refs,
    FilterConditions* filters,
    GenericFilters* gfilters,
    bint apply_filters
) except -1 nogil
