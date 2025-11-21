# cython: language_level=3
"""
Unified TSV writer with optional compression - C API declarations.
"""

from libc.stdio cimport FILE

cdef extern from "zlib.h":
    ctypedef void* gzFile

# Compression types
cdef enum CompressionType:
    COMPRESSION_NONE = 0
    COMPRESSION_PIGZ = 1
    COMPRESSION_ZLIB = 2

# Writer structure
ctypedef struct TSVWriter:
    FILE* file_handle
    gzFile gz_handle
    CompressionType compression_type
    int compression_level
    int compression_threads
    char* line_buffer
    size_t line_capacity
    size_t line_pos
    char* output_buffer
    size_t output_capacity
    size_t output_pos

# Writer lifecycle
cdef TSVWriter* tsv_writer_open(const char* path, int compression_level, int compression_threads) except NULL nogil
cdef int tsv_writer_flush(TSVWriter* writer) except -1 nogil
cdef int tsv_writer_close(TSVWriter* writer) except -1 nogil

# Row building
cdef int tsv_row_start(TSVWriter* writer) except -1 nogil
cdef int tsv_append_string(TSVWriter* writer, const char* s) except -1 nogil
cdef int tsv_append_int(TSVWriter* writer, long long val) except -1 nogil
cdef int tsv_append_float(TSVWriter* writer, double val, int precision) except -1 nogil
cdef int tsv_row_end(TSVWriter* writer) except -1 nogil

# Raw write (for headers)
cdef int tsv_write_raw(TSVWriter* writer, const char* data, size_t length) except -1 nogil
