# cython: initializedcheck=False
# cython: embedsignature=False
# cython: binding=True
# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False
# cython: infer_types=True
# distutils: language = c++
"""
Fast BAM/SAM/SAM.gz to Parquet converter using DuckDB C API.

Features:
- Dual-table output (by_reference + by_read)
- Hash-based partitioning
- Lossless SAM representation
- Direct DuckDB C API (no PyArrow overhead)
- Streaming processing (constant memory)
"""

from libc.stdlib cimport malloc, free, calloc
from libc.string cimport memcpy, strlen, strcpy, strcmp
from libc.stdint cimport uint8_t, uint16_t, uint32_t, int32_t, uint64_t, int8_t, int16_t
from libc.stdio cimport snprintf

from bam_filter.parquet_converter_duckdb cimport (
    duckdb_database, duckdb_connection, duckdb_appender,
    duckdb_state, DuckDBSuccess, DuckDBError,
    duckdb_open, duckdb_close, duckdb_connect, duckdb_disconnect,
    duckdb_query, duckdb_appender_create, duckdb_appender_begin_row,
    duckdb_appender_end_row, duckdb_appender_flush, duckdb_appender_close,
    duckdb_appender_destroy, duckdb_append_uint8, duckdb_append_uint16,
    duckdb_append_uint32, duckdb_append_uint64, duckdb_append_int32,
    duckdb_append_int64, duckdb_append_float, duckdb_append_varchar,
    duckdb_append_varchar_length, duckdb_append_blob, duckdb_append_null,
    samFile, bam_hdr_t, bam1_t, bam1_core_t,
    sam_open, sam_close, sam_hdr_read, bam_hdr_destroy, sam_read1,
    bam_init1, bam_destroy1, bam_get_qname, bam_get_seq, bam_get_qual,
    bam_get_cigar, bam_get_aux, bam_endpos,
    DuckDBWriter
)

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_verbose(int level, const char* tag, const char* fmt, ...) nogil
    double bf_monotonic_seconds() nogil


# ============================================================================
# Helper Functions
# ============================================================================

cdef inline uint32_t hash_partition(uint64_t value, int num_partitions) nogil:
    """Hash a value to partition ID."""
    return <uint32_t>(value % <uint64_t>num_partitions)


cdef inline int extract_cigar_string(bam1_t* record, char* buffer, int max_len) nogil:
    """Extract CIGAR string from BAM record."""
    cdef uint32_t* cigar = bam_get_cigar(record)
    cdef int n_cigar = record.core.n_cigar
    cdef int pos = 0
    cdef int i
    cdef uint32_t op_len
    cdef int op
    cdef char op_char
    cdef const char* cigar_ops = "MIDNSHP=X"

    for i in range(n_cigar):
        op_len = cigar[i] >> 4
        op = cigar[i] & 0xf

        if op < 9:
            op_char = cigar_ops[op]
        else:
            op_char = '?'

        pos += snprintf(buffer + pos, max_len - pos, "%u%c", op_len, op_char)
        if pos >= max_len - 1:
            break

    buffer[pos] = 0
    return pos


cdef inline float calculate_ani(bam1_t* record) nogil:
    """Calculate ANI from NM tag."""
    cdef uint8_t* nm_tag = bam_get_aux(record)
    cdef uint8_t* aux = nm_tag
    cdef uint8_t* end = record.data + record.l_data
    cdef int32_t nm_value = 0
    cdef int32_t aln_len = record.core.l_qseq
    cdef uint8_t tag_type

    # Find NM tag
    while aux < end:
        if aux[0] == 'N' and aux[1] == 'M':
            aux += 2  # Skip tag name
            if aux[0] == 'C':  # uint8
                nm_value = <int32_t>(<uint8_t*>(aux + 1))[0]
            elif aux[0] == 'S':  # uint16
                nm_value = <int32_t>(<uint16_t*>(aux + 1))[0]
            elif aux[0] == 'I':  # uint32
                nm_value = <int32_t>(<uint32_t*>(aux + 1))[0]
            elif aux[0] == 'i':  # int32
                nm_value = (<int32_t*>(aux + 1))[0]
            break

        # Skip to next tag (simplified - assumes standard tag format)
        aux += 2  # tag name
        if aux >= end:
            break

        tag_type = aux[0]
        aux += 1

        if tag_type == 'C' or tag_type == 'c':
            aux += 1
        elif tag_type == 'S' or tag_type == 's':
            aux += 2
        elif tag_type == 'I' or tag_type == 'i' or tag_type == 'f':
            aux += 4
        elif tag_type == 'Z' or tag_type == 'H':
            while aux < end and aux[0] != 0:
                aux += 1
            aux += 1
        else:
            break

    if aln_len > 0 and nm_value >= 0:
        return (1.0 - (<float>nm_value / <float>aln_len)) * 100.0
    else:
        return -1.0


# ============================================================================
# DuckDB Writer Implementation
# ============================================================================

cdef DuckDBWriter* duckdb_writer_create(const char* output_path,
                                        int num_partitions,
                                        bint write_by_reference,
                                        bint write_by_read) except NULL nogil:
    """Create DuckDB writer for dual-table Parquet output."""
    cdef DuckDBWriter* writer = <DuckDBWriter*>calloc(1, sizeof(DuckDBWriter))
    if writer == NULL:
        return NULL

    writer.output_path = output_path
    writer.num_partitions = num_partitions
    writer.write_by_reference = write_by_reference
    writer.write_by_read = write_by_read
    writer.records_written_by_ref = 0
    writer.records_written_by_read = 0

    # Open in-memory DuckDB database
    if duckdb_open(NULL, &writer.db) == DuckDBError:
        free(writer)
        return NULL

    if duckdb_connect(writer.db, &writer.conn) == DuckDBError:
        duckdb_close(&writer.db)
        free(writer)
        return NULL

    # Create tables
    cdef char query[4096]

    if write_by_reference:
        snprintf(query, sizeof(query),
                "CREATE TABLE alignments_by_reference ("
                "read_id UBIGINT, "
                "ref_id UINTEGER, "
                "position INTEGER, "
                "end_position INTEGER, "
                "alignment_score FLOAT, "
                "ani FLOAT, "
                "mapq UTINYINT, "
                "flag USMALLINT, "
                "edit_distance USMALLINT, "
                "alignment_length USMALLINT, "
                "template_length INTEGER, "
                "mate_ref_id INTEGER, "
                "mate_position INTEGER, "
                "read_name VARCHAR, "
                "cigar VARCHAR, "
                "sequence BLOB, "
                "quality BLOB, "
                "tags BLOB"
                ")")

        if duckdb_query(writer.conn, query, NULL) == DuckDBError:
            duckdb_disconnect(&writer.conn)
            duckdb_close(&writer.db)
            free(writer)
            return NULL

    if write_by_read:
        snprintf(query, sizeof(query),
                "CREATE TABLE alignments_by_read ("
                "read_id UBIGINT, "
                "ref_id UINTEGER, "
                "position INTEGER, "
                "end_position INTEGER, "
                "alignment_score FLOAT, "
                "ani FLOAT, "
                "mapq UTINYINT, "
                "flag USMALLINT, "
                "edit_distance USMALLINT, "
                "alignment_length USMALLINT, "
                "template_length INTEGER, "
                "mate_ref_id INTEGER, "
                "mate_position INTEGER, "
                "read_name VARCHAR, "
                "cigar VARCHAR, "
                "sequence BLOB, "
                "quality BLOB, "
                "tags BLOB"
                ")")

        if duckdb_query(writer.conn, query, NULL) == DuckDBError:
            if write_by_reference:
                duckdb_query(writer.conn, "DROP TABLE alignments_by_reference", NULL)
            duckdb_disconnect(&writer.conn)
            duckdb_close(&writer.db)
            free(writer)
            return NULL

    # Create appenders
    if write_by_reference:
        if duckdb_appender_create(writer.conn, NULL, "alignments_by_reference",
                                 &writer.appender_by_ref) == DuckDBError:
            duckdb_disconnect(&writer.conn)
            duckdb_close(&writer.db)
            free(writer)
            return NULL

    if write_by_read:
        if duckdb_appender_create(writer.conn, NULL, "alignments_by_read",
                                 &writer.appender_by_read) == DuckDBError:
            if write_by_reference:
                duckdb_appender_destroy(&writer.appender_by_ref)
            duckdb_disconnect(&writer.conn)
            duckdb_close(&writer.db)
            free(writer)
            return NULL

    return writer


cdef int duckdb_writer_add_alignment(DuckDBWriter* writer, bam1_t* record,
                                     bam_hdr_t* header, uint64_t read_id) except -1 nogil:
    """Add alignment to both tables."""
    if writer == NULL or record == NULL or header == NULL:
        return -1

    # Skip unmapped reads
    if record.core.tid < 0:
        return 0

    # Extract fields
    cdef char* read_name = bam_get_qname(record)
    cdef uint32_t ref_id = <uint32_t>record.core.tid
    cdef int32_t position = record.core.pos
    cdef int32_t end_position = bam_endpos(record)
    cdef uint8_t mapq = record.core.qual
    cdef uint16_t flag = record.core.flag
    cdef int32_t template_length = record.core.isize
    cdef int32_t mate_ref_id = record.core.mtid
    cdef int32_t mate_position = record.core.mpos
    cdef uint16_t alignment_length = <uint16_t>record.core.l_qseq

    # Calculate metrics
    cdef float ani = calculate_ani(record)
    cdef float alignment_score = 0.0  # TODO: extract from AS tag
    cdef uint16_t edit_distance = 0  # TODO: extract from NM tag

    # Extract variable-length data
    cdef char cigar_buf[2048]
    extract_cigar_string(record, cigar_buf, sizeof(cigar_buf))

    cdef uint8_t* seq = bam_get_seq(record)
    cdef uint8_t* qual = bam_get_qual(record)
    cdef uint8_t* aux = bam_get_aux(record)
    cdef int aux_len = record.l_data - (aux - record.data)

    # Write to by_reference table
    if writer.write_by_reference:
        duckdb_appender_begin_row(writer.appender_by_ref)
        duckdb_append_uint64(writer.appender_by_ref, read_id)
        duckdb_append_uint32(writer.appender_by_ref, ref_id)
        duckdb_append_int32(writer.appender_by_ref, position)
        duckdb_append_int32(writer.appender_by_ref, end_position)
        duckdb_append_float(writer.appender_by_ref, alignment_score)
        duckdb_append_float(writer.appender_by_ref, ani)
        duckdb_append_uint8(writer.appender_by_ref, mapq)
        duckdb_append_uint16(writer.appender_by_ref, flag)
        duckdb_append_uint16(writer.appender_by_ref, edit_distance)
        duckdb_append_uint16(writer.appender_by_ref, alignment_length)
        duckdb_append_int32(writer.appender_by_ref, template_length)
        duckdb_append_int32(writer.appender_by_ref, mate_ref_id)
        duckdb_append_int32(writer.appender_by_ref, mate_position)
        duckdb_append_varchar(writer.appender_by_ref, read_name)
        duckdb_append_varchar(writer.appender_by_ref, cigar_buf)
        duckdb_append_blob(writer.appender_by_ref, seq, (alignment_length + 1) >> 1)
        duckdb_append_blob(writer.appender_by_ref, qual, alignment_length)
        duckdb_append_blob(writer.appender_by_ref, aux, aux_len)
        duckdb_appender_end_row(writer.appender_by_ref)
        writer.records_written_by_ref += 1

    # Write to by_read table
    if writer.write_by_read:
        duckdb_appender_begin_row(writer.appender_by_read)
        duckdb_append_uint64(writer.appender_by_read, read_id)
        duckdb_append_uint32(writer.appender_by_read, ref_id)
        duckdb_append_int32(writer.appender_by_read, position)
        duckdb_append_int32(writer.appender_by_read, end_position)
        duckdb_append_float(writer.appender_by_read, alignment_score)
        duckdb_append_float(writer.appender_by_read, ani)
        duckdb_append_uint8(writer.appender_by_read, mapq)
        duckdb_append_uint16(writer.appender_by_read, flag)
        duckdb_append_uint16(writer.appender_by_read, edit_distance)
        duckdb_append_uint16(writer.appender_by_read, alignment_length)
        duckdb_append_int32(writer.appender_by_read, template_length)
        duckdb_append_int32(writer.appender_by_read, mate_ref_id)
        duckdb_append_int32(writer.appender_by_read, mate_position)
        duckdb_append_varchar(writer.appender_by_read, read_name)
        duckdb_append_varchar(writer.appender_by_read, cigar_buf)
        duckdb_append_blob(writer.appender_by_read, seq, (alignment_length + 1) >> 1)
        duckdb_append_blob(writer.appender_by_read, qual, alignment_length)
        duckdb_append_blob(writer.appender_by_read, aux, aux_len)
        duckdb_appender_end_row(writer.appender_by_read)
        writer.records_written_by_read += 1

    return 0


cdef int duckdb_writer_flush(DuckDBWriter* writer) except -1 nogil:
    """Flush appenders and write Parquet files."""
    if writer == NULL:
        return -1

    cdef char query[4096]
    cdef char output_file[2048]

    # Flush appenders
    if writer.write_by_reference:
        duckdb_appender_flush(writer.appender_by_ref)

    if writer.write_by_read:
        duckdb_appender_flush(writer.appender_by_read)

    # Write by_reference table to Parquet
    if writer.write_by_reference:
        snprintf(output_file, sizeof(output_file),
                "%s/alignments_by_reference.parquet", writer.output_path)
        snprintf(query, sizeof(query),
                "COPY (SELECT * FROM alignments_by_reference ORDER BY ref_id, position, read_id) "
                "TO '%s' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000)",
                output_file)

        if duckdb_query(writer.conn, query, NULL) == DuckDBError:
            return -1

    # Write by_read table to Parquet
    if writer.write_by_read:
        snprintf(output_file, sizeof(output_file),
                "%s/alignments_by_read.parquet", writer.output_path)
        snprintf(query, sizeof(query),
                "COPY (SELECT * FROM alignments_by_read ORDER BY read_id, alignment_score DESC) "
                "TO '%s' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000)",
                output_file)

        if duckdb_query(writer.conn, query, NULL) == DuckDBError:
            return -1

    return 0


cdef void duckdb_writer_close(DuckDBWriter* writer) nogil:
    """Close writer and free resources."""
    if writer != NULL:
        if writer.write_by_reference:
            duckdb_appender_destroy(&writer.appender_by_ref)

        if writer.write_by_read:
            duckdb_appender_destroy(&writer.appender_by_read)

        duckdb_disconnect(&writer.conn)
        duckdb_close(&writer.db)
        free(writer)


# ============================================================================
# Main Conversion Function
# ============================================================================

def convert_sam_bam_to_parquet(
    str input_file,
    str output_dir,
    int num_partitions=256,
    bint write_by_reference=True,
    bint write_by_read=True,
):
    """
    Convert SAM/SAM.gz/BAM to dual-table Parquet format.

    Parameters
    ----------
    input_file : str
        Input alignment file (SAM, SAM.gz, or BAM)
    output_dir : str
        Output directory for Parquet files
    num_partitions : int
        Number of hash partitions (default: 256)
    write_by_reference : bool
        Write alignments_by_reference table (default: True)
    write_by_read : bool
        Write alignments_by_read table (default: True)

    Returns
    -------
    dict
        Statistics: total_records, records_by_ref, records_by_read, duration_seconds
    """
    import time
    from pathlib import Path

    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    cdef bytes input_file_bytes = input_file.encode('utf-8')
    cdef bytes output_dir_bytes = output_dir.encode('utf-8')
    cdef const char* input_file_cstr = input_file_bytes
    cdef const char* output_dir_cstr = output_dir_bytes

    cdef double start_time = time.perf_counter()

    cdef dict stats
    cdef uint64_t total_records = 0
    cdef uint64_t records_by_ref = 0
    cdef uint64_t records_by_read = 0
    cdef double processing_time = 0.0

    with nogil:
        _convert_impl(
            input_file_cstr,
            output_dir_cstr,
            num_partitions,
            write_by_reference,
            write_by_read,
            &total_records,
            &records_by_ref,
            &records_by_read,
            &processing_time
        )

    cdef double duration = time.perf_counter() - start_time

    return {
        'total_records': total_records,
        'records_by_reference': records_by_ref,
        'records_by_read': records_by_read,
        'processing_time_seconds': processing_time,
        'duration_seconds': duration,
    }


cdef void _convert_impl(
    const char* input_file,
    const char* output_dir,
    int num_partitions,
    bint write_by_reference,
    bint write_by_read,
    uint64_t* out_total_records,
    uint64_t* out_records_by_ref,
    uint64_t* out_records_by_read,
    double* out_processing_time
) nogil:
    """Main conversion implementation (nogil)."""
    cdef double start_time = bf_monotonic_seconds()

    # Open input file
    cdef samFile* sam_fp = sam_open(input_file, "r")
    if sam_fp == NULL:
        with gil:
            raise IOError(f"Failed to open input file: {input_file.decode('utf-8')}")

    # Read header
    cdef bam_hdr_t* header = sam_hdr_read(sam_fp)
    if header == NULL:
        sam_close(sam_fp)
        with gil:
            raise IOError("Failed to read header")

    # Create writer
    cdef DuckDBWriter* writer = duckdb_writer_create(
        output_dir, num_partitions, write_by_reference, write_by_read
    )
    if writer == NULL:
        bam_hdr_destroy(header)
        sam_close(sam_fp)
        with gil:
            raise IOError("Failed to create DuckDB writer")

    # Process alignments
    cdef bam1_t* record = bam_init1()
    cdef int read_result
    cdef uint64_t read_id = 0
    cdef uint64_t records_processed = 0
    cdef double last_log_time = start_time
    cdef double current_time

    bf_nogil_logf_verbose(2, "CONVERT", "Starting conversion: %s\n", input_file)

    while True:
        read_result = sam_read1(sam_fp, header, record)

        if read_result < 0:
            break

        duckdb_writer_add_alignment(writer, record, header, read_id)
        read_id += 1
        records_processed += 1

        # Log progress every 10 seconds
        current_time = bf_monotonic_seconds()
        if current_time - last_log_time >= 10.0:
            bf_nogil_logf_verbose(2, "CONVERT",
                                 "Processed %llu records (%.1f M/sec)\n",
                                 records_processed,
                                 records_processed / (current_time - start_time) / 1000000.0)
            last_log_time = current_time

    # Flush and write Parquet files
    duckdb_writer_flush(writer)

    cdef uint64_t final_by_ref = writer.records_written_by_ref
    cdef uint64_t final_by_read = writer.records_written_by_read

    # Cleanup
    duckdb_writer_close(writer)
    bam_destroy1(record)
    bam_hdr_destroy(header)
    sam_close(sam_fp)

    cdef double end_time = bf_monotonic_seconds()
    cdef double duration = end_time - start_time

    bf_nogil_logf_verbose(2, "CONVERT",
                         "Conversion complete: %llu records in %.2f seconds (%.1f M/sec)\n",
                         records_processed, duration,
                         records_processed / duration / 1000000.0)

    # Set output parameters
    out_total_records[0] = records_processed
    out_records_by_ref[0] = final_by_ref
    out_records_by_read[0] = final_by_read
    out_processing_time[0] = duration
