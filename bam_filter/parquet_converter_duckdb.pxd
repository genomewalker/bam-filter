# cython: language_level=3
"""
DuckDB-based Parquet converter - C API declarations.
Uses DuckDB C++ API to write Parquet files without PyArrow.
"""

from libc.stdint cimport uint8_t, uint16_t, uint32_t, int32_t, uint64_t, int64_t

# DuckDB C API declarations
cdef extern from "duckdb.h":
    ctypedef struct duckdb_database:
        pass

    ctypedef struct duckdb_connection:
        pass

    ctypedef struct duckdb_appender:
        pass

    ctypedef enum duckdb_type:
        DUCKDB_TYPE_BOOLEAN = 0
        DUCKDB_TYPE_TINYINT = 1
        DUCKDB_TYPE_SMALLINT = 2
        DUCKDB_TYPE_INTEGER = 3
        DUCKDB_TYPE_BIGINT = 4
        DUCKDB_TYPE_UTINYINT = 5
        DUCKDB_TYPE_USMALLINT = 6
        DUCKDB_TYPE_UINTEGER = 7
        DUCKDB_TYPE_UBIGINT = 8
        DUCKDB_TYPE_FLOAT = 9
        DUCKDB_TYPE_DOUBLE = 10
        DUCKDB_TYPE_VARCHAR = 13
        DUCKDB_TYPE_BLOB = 14

    ctypedef enum duckdb_state:
        DuckDBSuccess = 0
        DuckDBError = 1

    # Database operations
    duckdb_state duckdb_open(const char* path, duckdb_database* out_database) nogil
    void duckdb_close(duckdb_database* database) nogil
    duckdb_state duckdb_connect(duckdb_database database, duckdb_connection* out_connection) nogil
    void duckdb_disconnect(duckdb_connection* connection) nogil

    # Query operations
    duckdb_state duckdb_query(duckdb_connection connection, const char* query, void* out_result) nogil

    # Appender operations
    duckdb_state duckdb_appender_create(duckdb_connection connection,
                                        const char* schema,
                                        const char* table,
                                        duckdb_appender* out_appender) nogil
    duckdb_state duckdb_appender_begin_row(duckdb_appender appender) nogil
    duckdb_state duckdb_appender_end_row(duckdb_appender appender) nogil
    duckdb_state duckdb_appender_flush(duckdb_appender appender) nogil
    duckdb_state duckdb_appender_close(duckdb_appender appender) nogil
    duckdb_state duckdb_appender_destroy(duckdb_appender* appender) nogil

    # Append data
    duckdb_state duckdb_append_uint8(duckdb_appender appender, uint8_t value) nogil
    duckdb_state duckdb_append_uint16(duckdb_appender appender, uint16_t value) nogil
    duckdb_state duckdb_append_uint32(duckdb_appender appender, uint32_t value) nogil
    duckdb_state duckdb_append_uint64(duckdb_appender appender, uint64_t value) nogil
    duckdb_state duckdb_append_int32(duckdb_appender appender, int32_t value) nogil
    duckdb_state duckdb_append_int64(duckdb_appender appender, int64_t value) nogil
    duckdb_state duckdb_append_float(duckdb_appender appender, float value) nogil
    duckdb_state duckdb_append_double(duckdb_appender appender, double value) nogil
    duckdb_state duckdb_append_varchar(duckdb_appender appender, const char* value) nogil
    duckdb_state duckdb_append_varchar_length(duckdb_appender appender, const char* value, uint64_t length) nogil
    duckdb_state duckdb_append_blob(duckdb_appender appender, const void* data, uint64_t length) nogil
    duckdb_state duckdb_append_null(duckdb_appender appender) nogil

# HTSlib declarations
cdef extern from "htslib/sam.h":
    ctypedef struct samFile:
        pass

    ctypedef struct bam_hdr_t:
        int32_t n_targets
        char** target_name
        uint32_t* target_len

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

    char* bam_get_qname(bam1_t* b) nogil
    uint8_t* bam_get_seq(bam1_t* b) nogil
    uint8_t* bam_get_qual(bam1_t* b) nogil
    uint32_t* bam_get_cigar(bam1_t* b) nogil
    uint8_t* bam_get_aux(bam1_t* b) nogil

    int bam_endpos(bam1_t* b) nogil

# Converter structures
ctypedef struct DuckDBWriter:
    duckdb_database db
    duckdb_connection conn
    duckdb_appender appender_by_ref
    duckdb_appender appender_by_read
    const char* output_path
    int num_partitions
    uint64_t records_written_by_ref
    uint64_t records_written_by_read
    bint write_by_reference
    bint write_by_read

# Converter API
cdef DuckDBWriter* duckdb_writer_create(const char* output_path,
                                        int num_partitions,
                                        bint write_by_reference,
                                        bint write_by_read) except NULL nogil
cdef int duckdb_writer_add_alignment(DuckDBWriter* writer, bam1_t* record,
                                     bam_hdr_t* header, uint64_t read_id) except -1 nogil
cdef int duckdb_writer_flush(DuckDBWriter* writer) except -1 nogil
cdef void duckdb_writer_close(DuckDBWriter* writer) nogil
