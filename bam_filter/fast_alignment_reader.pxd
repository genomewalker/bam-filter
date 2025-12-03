# cython: language_level=3
"""
Fast alignment reader - C API declarations.
Unified reader for SAM, SAM.gz, and BAM formats.
"""

from libc.stdint cimport uint8_t, uint16_t, uint32_t, int32_t, uint64_t
from libc.stdio cimport FILE

# HTSlib declarations
cdef extern from "htslib/sam.h":
    ctypedef struct samFile:
        pass

    ctypedef struct bam_hdr_t:
        int32_t n_targets
        char** target_name
        uint32_t* target_len
        char* text
        uint32_t l_text

    ctypedef struct bam1_core_t:
        int32_t pos
        int32_t tid
        uint16_t bin
        uint8_t qual
        uint8_t l_qname
        uint16_t flag
        uint16_t n_cigar
        int32_t l_qseq
        int32_t mtid
        int32_t mpos
        int32_t isize

    ctypedef struct bam1_t:
        bam1_core_t core
        int l_data
        int m_data
        uint8_t* data

    samFile* sam_open(const char* fn, const char* mode) nogil
    int sam_close(samFile* fp) nogil
    bam_hdr_t* sam_hdr_read(samFile* fp) nogil
    void bam_hdr_destroy(bam_hdr_t* h) nogil
    int sam_read1(samFile* fp, bam_hdr_t* h, bam1_t* b) nogil
    bam1_t* bam_init1() nogil
    void bam_destroy1(bam1_t* b) nogil

    uint8_t* bam_get_seq(bam1_t* b) nogil
    uint8_t* bam_get_qual(bam1_t* b) nogil
    char* bam_get_qname(bam1_t* b) nogil
    uint32_t* bam_get_cigar(bam1_t* b) nogil
    uint8_t* bam_get_aux(bam1_t* b) nogil

# Format detection
cdef enum AlignmentFormat:
    FORMAT_UNKNOWN = 0
    FORMAT_SAM = 1
    FORMAT_SAM_GZ = 2
    FORMAT_BAM = 3

# Alignment reader structure
ctypedef struct AlignmentReader:
    AlignmentFormat format
    samFile* hts_file
    bam_hdr_t* header
    bam1_t* current_record
    uint64_t records_read
    bint is_open

# Reader API
cdef AlignmentReader* alignment_reader_open(const char* path) except NULL nogil
cdef int alignment_reader_read(AlignmentReader* reader, bam1_t** record) except -1 nogil
cdef int alignment_reader_close(AlignmentReader* reader) except -1 nogil
cdef AlignmentFormat detect_format(const char* path) nogil

# Helper functions for field extraction
cdef int extract_read_name(bam1_t* record, char* buffer, int max_len) nogil
cdef int extract_sequence(bam1_t* record, char* buffer, int max_len) nogil
cdef int extract_quality(bam1_t* record, char* buffer, int max_len) nogil
cdef int extract_cigar_string(bam1_t* record, char* buffer, int max_len) nogil
cdef uint8_t* find_tag(bam1_t* record, const char* tag) nogil
cdef int32_t get_tag_int(bam1_t* record, const char* tag, int32_t default_val) nogil
cdef float get_tag_float(bam1_t* record, const char* tag, float default_val) nogil
