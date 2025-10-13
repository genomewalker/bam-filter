# cython: language_level=3

from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t

# Import types from processor.pxd so MemoryPool/ReferenceMapping and bam/htslib types
# are available for declarations in this pxd.
from bam_filter.processor cimport MemoryPool, bam1_t, sam_hdr_t, samFile
from bam_filter.processor_mapping cimport ReferenceMapping, create_reference_mapping, destroy_reference_mapping

cdef struct CompactAlignment:
    uint32_t original_position
    uint32_t pool_index
    uint32_t reference_id

cdef struct LookupTable:
    CompactAlignment* alignments
    uint64_t* ref_starts
    uint32_t* ref_counts
    uint64_t total_count
    uint32_t num_refs

cdef struct WriteBatch:
    bam1_t** records
    uint64_t* pool_indices
    uint32_t count
    uint32_t capacity
    uint32_t reference_id

# Expose functions implemented in processor_bam_writer.pyx
cdef LookupTable* create_lookup_table(MemoryPool* pool, ReferenceMapping* mapping) except NULL nogil
cdef void destroy_lookup_table(LookupTable* table) noexcept nogil

cdef WriteBatch* create_write_batch(uint32_t capacity) except NULL nogil
cdef void destroy_write_batch(WriteBatch* batch) noexcept nogil

# This function performs higher-level I/O and uses Python/C-API safe
# operations; do not declare it nogil so callers will hold the GIL.
cdef int write_filtered_bam(MemoryPool* pool,
                             const char* input_bam_path,
                             const char* output_bam_path,
                             sam_hdr_t* header,
                             ReferenceMapping* mapping,
                             int num_threads) except -1

# Minimal HTSlib cimports (so other modules can use these types)
cdef int copy_bam_record(bam1_t* src, bam1_t* dst, int extra_bytes) except -1 nogil
cdef int write_batch_to_bam(samFile* out_bam, sam_hdr_t* header, WriteBatch* batch) except -1 nogil

cdef int write_bam_with_filtered_header(MemoryPool* pool,
                                       const char* input_bam_path,
                                       const char* output_bam_path,
                                       sam_hdr_t* original_header,
                                       ReferenceMapping* existing_mapping,
                                       int bam_write_threads) except -1
