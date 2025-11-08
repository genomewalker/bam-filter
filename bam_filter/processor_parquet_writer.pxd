# cython: language_level=3
# -*- coding: utf-8 -*-

"""BAM to Parquet converter C declarations."""

from libc.stdint cimport uint32_t, uint16_t, uint8_t, int32_t
from bam_filter.processor cimport bam1_t, sam_hdr_t, AlignmentScoringConfig


cdef struct ParquetAlignment:
    uint32_t read_id
    uint32_t ref_id
    int32_t position
    int32_t end_position
    uint8_t mapq
    uint16_t flag
    float ani
    float alignment_score
    float pmd_score
    uint32_t num_mismatches
    uint32_t alignment_length
    int32_t template_length
    int32_t mate_ref_id
    int32_t mate_position
    char* read_name
    char* cigar
    char* sequence
    uint8_t* quality
    uint32_t quality_length
    uint8_t* tags
    uint32_t tags_length


cdef char* extract_cigar_string(bam1_t* b) nogil
cdef char* extract_sequence_string(bam1_t* b) nogil
cdef int populate_parquet_alignment(bam1_t* b, sam_hdr_t* header,
                                    AlignmentScoringConfig* config,
                                    uint32_t read_id, ParquetAlignment* aln,
                                    bint include_read_name,
                                    bint include_sequence) except -1 nogil
cdef void free_parquet_alignment(ParquetAlignment* aln) noexcept nogil
cpdef object create_references_table(str bam_path, int num_partitions)
