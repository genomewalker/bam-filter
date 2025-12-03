# cython: language_level=3
"""
Fast dual-table Parquet writer - C API declarations.
Writes alignments partitioned by reference and by read for optimal query performance.
"""

from libc.stdint cimport uint8_t, uint16_t, uint32_t, int32_t, uint64_t, int64_t

# Alignment record structure (in-memory representation)
ctypedef struct AlignmentRecord:
    # IDs
    uint64_t read_id
    uint32_t ref_id

    # Position
    int32_t position
    int32_t end_position

    # Quality
    uint8_t mapq
    uint16_t flag

    # Scores and metrics
    float alignment_score
    float ani
    uint16_t edit_distance
    uint16_t alignment_length
    uint16_t reference_span

    # Calculated metrics (optional)
    float pmd_score
    float gc_content
    float dust_score

    # Paired-end
    int32_t template_length
    int32_t mate_ref_id
    int32_t mate_position

    # Variable-length data (pointers to external buffers)
    char* read_name
    char* cigar
    char* sequence
    char* quality
    uint8_t* tags
    int tags_len

    # Lengths
    uint16_t read_name_len
    uint16_t cigar_len
    uint16_t sequence_len
    uint16_t quality_len


# Batch writer for accumulating records before writing
ctypedef struct ParquetBatch:
    AlignmentRecord* records
    int64_t capacity
    int64_t count
    int partition_id


# Parquet writer configuration
ctypedef struct ParquetWriterConfig:
    const char* output_path
    int num_partitions
    int batch_size
    const char* compression
    int compression_level
    bint write_by_reference
    bint write_by_read
    bint include_read_names
    bint include_sequences
    int num_threads


# Parquet writer handle
ctypedef struct ParquetWriter:
    ParquetWriterConfig config

    # Statistics
    uint64_t total_records_written
    uint64_t total_bytes_written
    uint64_t records_by_reference
    uint64_t records_by_read

    # Per-partition batches (for by_reference table)
    ParquetBatch** ref_batches

    # Per-partition batches (for by_read table)
    ParquetBatch** read_batches

    # Reference metadata
    int num_references
    char** reference_names
    uint32_t* reference_lengths


# Writer API
cdef ParquetWriter* parquet_writer_create(ParquetWriterConfig* config) except NULL nogil
cdef int parquet_writer_add_record(ParquetWriter* writer, AlignmentRecord* record) except -1 nogil
cdef int parquet_writer_flush(ParquetWriter* writer) except -1 nogil
cdef int parquet_writer_close(ParquetWriter* writer) except -1 nogil

# Batch management
cdef ParquetBatch* create_batch(int capacity, int partition_id) except NULL nogil
cdef int batch_add_record(ParquetBatch* batch, AlignmentRecord* record) except -1 nogil
cdef int batch_write_to_parquet(ParquetBatch* batch, const char* output_file,
                                 ParquetWriterConfig* config, bint by_read) except -1 nogil
cdef void batch_clear(ParquetBatch* batch) nogil
cdef void batch_destroy(ParquetBatch* batch) nogil

# Helper functions
cdef int calculate_partition_id(uint64_t id_value, int num_partitions) nogil
cdef float calculate_ani(const char* cigar, int edit_distance, int alignment_length) nogil
cdef float calculate_gc_content(const char* sequence, int seq_len) nogil
