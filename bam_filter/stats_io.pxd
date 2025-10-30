# cython: language_level=3
from libc.stdio cimport FILE
from bam_filter.stats cimport RefStats, FilterConditions
from bam_filter.processor_types cimport sam_hdr_t

# zlib gzFile is declared directly to avoid depending on a non-portable
# Cython-provided zlib.pxd. Expose the minimal API we need for callers.
cdef extern from "zlib.h":
	ctypedef void* gzFile
	gzFile gzopen(const char* path, const char* mode) nogil
	int gzclose(gzFile file) nogil
	int gzprintf(gzFile file, const char* format, ...) nogil

# Declare the IO functions exported by stats_io. They are nogil-capable
# and return int; provide an explicit exception value to avoid forcing
# exception checks on cimporting callers (performance hint).
cdef int write_stats_to_file(FILE* fp, RefStats* global_ref_stats, sam_hdr_t* header, int n_refs, FilterConditions* filters, bint apply_filters) except -1 nogil

cdef int write_stats_to_gzfile(gzFile gzfp, RefStats* global_ref_stats, sam_hdr_t* header, int n_refs, FilterConditions* filters, bint apply_filters) except -1 nogil

cdef int write_output_files_complete(const char* output_c, const char* filtered_output_c, RefStats* global_ref_stats, sam_hdr_t* header, int n_refs, FilterConditions* filters) except -1 nogil
