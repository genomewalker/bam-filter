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
Streaming Parquet converter with parallel processing and batched writes.

Key improvements:
- Parallel SAM/BAM reading using OpenMP
- Individual tag columns (AS, NM, MD) for better compression
- Batched streaming writes every 100K records
- Hash-partitioned output files
- Minimal memory footprint
"""

from libc.stdlib cimport malloc, free, calloc
from libc.string cimport memcpy, strlen, strcpy, strcmp, memset
from libc.stdint cimport uint8_t, uint16_t, uint32_t, int32_t, uint64_t, int8_t, int16_t
from libc.stdio cimport snprintf, FILE, fopen, fclose, fprintf
from cython.parallel cimport prange, parallel
from libc.math cimport isnan

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
    bam_get_cigar, bam_get_aux, bam_endpos
)

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_verbose(int level, const char* tag, const char* fmt, ...) nogil
    double bf_monotonic_seconds() nogil

cdef extern from "omp.h":
    int omp_get_thread_num() nogil


# ============================================================================
# Improved Tag Extraction (individual columns for common tags)
# ============================================================================

cdef inline int32_t extract_tag_int(bam1_t* record, const char* tag_name) nogil:
    """Extract integer tag (AS, NM, etc.)."""
    cdef uint8_t* aux = bam_get_aux(record)
    cdef uint8_t* end = record.data + record.l_data
    cdef uint8_t tag_type

    while aux < end - 2:
        if aux[0] == tag_name[0] and aux[1] == tag_name[1]:
            aux += 2  # Skip tag name
            if aux[0] == b'C':  # uint8
                return <int32_t>(<uint8_t*>(aux + 1))[0]
            elif aux[0] == b'c':  # int8
                return <int32_t>(<int8_t*>(aux + 1))[0]
            elif aux[0] == b'S':  # uint16
                return <int32_t>(<uint16_t*>(aux + 1))[0]
            elif aux[0] == b's':  # int16
                return <int32_t>(<int16_t*>(aux + 1))[0]
            elif aux[0] == b'I':  # uint32
                return <int32_t>(<uint32_t*>(aux + 1))[0]
            elif aux[0] == b'i':  # int32
                return (<int32_t*>(aux + 1))[0]
            return -1

        # Skip to next tag
        aux += 2
        if aux >= end:
            break

        tag_type = aux[0]
        aux += 1

        if tag_type == b'C' or tag_type == b'c':
            aux += 1
        elif tag_type == b'S' or tag_type == b's':
            aux += 2
        elif tag_type == b'I' or tag_type == b'i' or tag_type == b'f':
            aux += 4
        elif tag_type == b'Z' or tag_type == b'H':
            while aux < end and aux[0] != 0:
                aux += 1
            aux += 1
        else:
            break

    return -1


cdef inline int extract_tag_string(bam1_t* record, const char* tag_name, char* buffer, int max_len) nogil:
    """Extract string tag (MD, etc.)."""
    cdef uint8_t* aux = bam_get_aux(record)
    cdef uint8_t* end = record.data + record.l_data
    cdef int len_written = 0
    cdef uint8_t tag_type

    buffer[0] = 0

    while aux < end - 2:
        if aux[0] == tag_name[0] and aux[1] == tag_name[1]:
            aux += 2  # Skip tag name
            if aux[0] == b'Z':  # String
                aux += 1  # Skip type
                while aux < end and aux[0] != 0 and len_written < max_len - 1:
                    buffer[len_written] = <char>aux[0]
                    len_written += 1
                    aux += 1
                buffer[len_written] = 0
                return len_written
            return 0

        # Skip to next tag
        aux += 2
        if aux >= end:
            break

        tag_type = aux[0]
        aux += 1

        if tag_type == b'C' or tag_type == b'c':
            aux += 1
        elif tag_type == b'S' or tag_type == b's':
            aux += 2
        elif tag_type == b'I' or tag_type == b'i' or tag_type == b'f':
            aux += 4
        elif tag_type == b'Z' or tag_type == b'H':
            while aux < end and aux[0] != 0:
                aux += 1
            aux += 1
        else:
            break

    return 0


cdef inline float extract_tag_float(bam1_t* record, const char* tag_name) nogil:
    """Extract float tag (PMD, ZS, ZP, etc.) - optimized manual parsing."""
    cdef uint8_t* aux = bam_get_aux(record)
    cdef uint8_t* end = record.data + record.l_data
    cdef uint8_t tag_type

    if aux == NULL or aux >= end - 2:
        return -1.0

    while aux < end - 2:
        if aux[0] == tag_name[0] and aux[1] == tag_name[1]:
            aux += 2  # Skip tag name
            tag_type = aux[0]
            aux += 1

            if tag_type == b'f':  # Float (4 bytes)
                return (<float*>aux)[0]
            elif tag_type == b'd':  # Double (8 bytes)
                return <float>((<double*>aux)[0])
            elif tag_type == b'i':  # int32
                return <float>((<int32_t*>aux)[0])
            elif tag_type == b'I':  # uint32
                return <float>((<uint32_t*>aux)[0])
            elif tag_type == b'c':  # int8
                return <float>((<int8_t*>aux)[0])
            elif tag_type == b'C':  # uint8
                return <float>((<uint8_t*>aux)[0])
            elif tag_type == b's':  # int16
                return <float>((<int16_t*>aux)[0])
            elif tag_type == b'S':  # uint16
                return <float>((<uint16_t*>aux)[0])
            return -1.0

        # Skip to next tag
        aux += 2
        if aux >= end:
            break

        tag_type = aux[0]
        aux += 1

        if tag_type == b'C' or tag_type == b'c':
            aux += 1
        elif tag_type == b'S' or tag_type == b's':
            aux += 2
        elif tag_type == b'I' or tag_type == b'i' or tag_type == b'f':
            aux += 4
        elif tag_type == b'd':
            aux += 8
        elif tag_type == b'Z' or tag_type == b'H':
            while aux < end and aux[0] != 0:
                aux += 1
            aux += 1  # Skip null terminator
        else:
            break

    return -1.0


# ============================================================================
# Partition Writer (writes batches to partition-specific DuckDB files)
# ============================================================================

ctypedef struct PartitionWriter:
    duckdb_database db
    duckdb_connection conn
    duckdb_appender appender
    char filepath[512]
    char tablename[128]
    uint64_t records_written
    uint32_t partition_id
    bint is_open


cdef PartitionWriter* partition_writer_create(const char* output_dir,
                                              uint32_t partition_id,
                                              const char* table_prefix,
                                              bint is_reference_table,
                                              int thread_id) except NULL nogil:
    """Create a partition-specific writer (per thread)."""
    cdef PartitionWriter* pw = <PartitionWriter*>calloc(1, sizeof(PartitionWriter))
    if pw == NULL:
        return NULL

    pw.partition_id = partition_id
    pw.records_written = 0
    pw.is_open = False

    # Create partition-specific file path with thread ID
    if is_reference_table:
        snprintf(pw.filepath, sizeof(pw.filepath),
                "%s/by_reference_p%04d_t%02d.duckdb", output_dir, partition_id, thread_id)
        snprintf(pw.tablename, sizeof(pw.tablename), "alignments")
    else:
        snprintf(pw.filepath, sizeof(pw.filepath),
                "%s/by_read_p%04d_t%02d.duckdb", output_dir, partition_id, thread_id)
        snprintf(pw.tablename, sizeof(pw.tablename), "alignments")

    # Open DuckDB file
    if duckdb_open(pw.filepath, &pw.db) == DuckDBError:
        free(pw)
        return NULL

    if duckdb_connect(pw.db, &pw.conn) == DuckDBError:
        duckdb_close(&pw.db)
        free(pw)
        return NULL

    # Create table with comprehensive pipeline state columns
    cdef char query[4096]
    snprintf(query, sizeof(query),
            "CREATE TABLE %s ("
            # Core identifiers
            "read_id UBIGINT, "
            "read_name VARCHAR, "
            "ref_id UINTEGER, "
            "ref_name VARCHAR, "

            # Alignment coordinates
            "position INTEGER, "
            "end_position INTEGER, "
            "mapq UTINYINT, "
            "flag USMALLINT, "
            "alignment_length USMALLINT, "

            # Paired-end info
            "template_length INTEGER, "
            "mate_ref_id INTEGER, "
            "mate_position INTEGER, "

            # Alignment quality (individual tag columns)
            "alignment_score INTEGER, "   # AS tag - primary aligner score
            "xs_score INTEGER, "          # XS tag - suboptimal alignment score
            "edit_distance USMALLINT, "   # NM tag - edit distance
            "num_mismatches UTINYINT, "   # XM tag - number of mismatches
            "num_gap_opens UTINYINT, "    # XO tag - gap opens
            "num_gap_extensions UTINYINT, " # XG tag - gap extensions
            "md_string VARCHAR, "         # MD tag - mismatch positions
            "ani FLOAT, "                 # Calculated: (1 - NM/length) * 100

            # Ancient DNA damage
            "pmd_score FLOAT, "           # PMD tag - ancient DNA damage probability

            # bam-filter pipeline state
            "zs_score FLOAT, "            # ZS tag - precomputed score (filter/reassign)
            "zp_posterior FLOAT, "        # ZP tag - EM posterior probability
            "lca_taxid INTEGER, "         # LCA-assigned taxonomic ID
            "reassigned_ref_id INTEGER, " # Reassigned reference ID (or -1 if unchanged)
            "filter_passed BOOLEAN, "     # Passed all filters (true) or filtered out (false)

            # Read group and metadata
            "read_group VARCHAR, "        # RG tag - read group identifier

            # Sequence data
            "cigar VARCHAR, "
            "sequence BLOB, "
            "quality BLOB"
            ")", pw.tablename)

    if duckdb_query(pw.conn, query, NULL) == DuckDBError:
        duckdb_disconnect(&pw.conn)
        duckdb_close(&pw.db)
        free(pw)
        return NULL

    # Create appender
    if duckdb_appender_create(pw.conn, NULL, pw.tablename, &pw.appender) == DuckDBError:
        duckdb_disconnect(&pw.conn)
        duckdb_close(&pw.db)
        free(pw)
        return NULL

    pw.is_open = True
    return pw


cdef int partition_writer_add_record(PartitionWriter* pw, bam1_t* record, bam_hdr_t* header,
                                     uint64_t read_id) except -1 nogil:
    """Add record to partition writer with comprehensive tag extraction."""
    if pw == NULL or not pw.is_open:
        bf_nogil_logf_verbose(1, "STREAM", "ERROR: Writer is NULL or not open\n")
        return -1

    # Debug: Log first few records
    if pw.records_written < 5:
        bf_nogil_logf_verbose(2, "STREAM", "Adding record %llu to partition %u (records_written=%llu)\n",
                            read_id, pw.partition_id, pw.records_written)

    # Extract basic fields
    cdef uint32_t ref_id = <uint32_t>record.core.tid
    cdef int32_t position = record.core.pos
    cdef int32_t end_position = bam_endpos(record)
    cdef uint8_t mapq = record.core.qual
    cdef uint16_t flag = record.core.flag
    cdef uint16_t alignment_length = <uint16_t>record.core.l_qseq
    cdef int32_t template_length = record.core.isize
    cdef int32_t mate_ref_id = record.core.mtid
    cdef int32_t mate_position = record.core.mpos

    # Extract read and reference names
    cdef char* read_name = bam_get_qname(record)
    cdef const char* ref_name

    # Safe reference name extraction with bounds checking
    if ref_id >= 0 and ref_id < header.n_targets and header.target_name != NULL:
        ref_name = header.target_name[ref_id]
        if ref_name == NULL:
            ref_name = "*"
    else:
        ref_name = "*"

    # Extract alignment quality tags
    cdef int32_t alignment_score = extract_tag_int(record, "AS")
    cdef int32_t xs_score = extract_tag_int(record, "XS")
    cdef int32_t edit_distance = extract_tag_int(record, "NM")
    cdef int32_t num_mismatches = extract_tag_int(record, "XM")
    cdef int32_t num_gap_opens = extract_tag_int(record, "XO")
    cdef int32_t num_gap_extensions = extract_tag_int(record, "XG")

    cdef char md_string[1024]
    cdef char read_group[256]
    extract_tag_string(record, "MD", md_string, sizeof(md_string))
    extract_tag_string(record, "RG", read_group, sizeof(read_group))

    # Extract pipeline state tags (float)
    cdef float pmd_score = extract_tag_float(record, "PM")  # or "PMD"
    if pmd_score < 0:
        pmd_score = extract_tag_float(record, "PMD")  # Try alternate tag name
    cdef float zs_score = extract_tag_float(record, "ZS")
    cdef float zp_posterior = extract_tag_float(record, "ZP")

    # Extract LCA taxid and reassignment info
    cdef int32_t lca_taxid = extract_tag_int(record, "ZT")  # LCA taxonomic ID
    cdef int32_t reassigned_ref_id = extract_tag_int(record, "ZR")  # Reassigned ref ID

    # Filter status (default to true if not present - unfiltered data)
    cdef int32_t filter_flag = extract_tag_int(record, "ZF")  # 1=passed, 0=filtered
    cdef bint filter_passed = (filter_flag != 0) if filter_flag >= 0 else True

    # Calculate ANI
    cdef float ani = -1.0
    if edit_distance >= 0 and alignment_length > 0:
        ani = (1.0 - (<float>edit_distance / <float>alignment_length)) * 100.0

    cdef char cigar_buf[2048]
    cdef uint32_t* cigar = bam_get_cigar(record)
    cdef int n_cigar = record.core.n_cigar
    cdef int pos = 0
    cdef const char* cigar_ops = "MIDNSHP=X"
    cdef int i
    cdef uint32_t op_len
    cdef int op
    cdef char op_char

    for i in range(n_cigar):
        op_len = cigar[i] >> 4
        op = cigar[i] & 0xf
        op_char = cigar_ops[op] if op < 9 else '?'
        pos += snprintf(cigar_buf + pos, sizeof(cigar_buf) - pos, "%u%c", op_len, op_char)
        if pos >= sizeof(cigar_buf) - 1:
            break
    cigar_buf[pos] = 0

    cdef uint8_t* seq = bam_get_seq(record)
    cdef uint8_t* qual = bam_get_qual(record)

    # Append row - must match schema order exactly
    duckdb_appender_begin_row(pw.appender)

    # Core identifiers
    duckdb_append_uint64(pw.appender, read_id)
    duckdb_append_varchar(pw.appender, read_name)
    duckdb_append_uint32(pw.appender, ref_id)
    duckdb_append_varchar(pw.appender, ref_name)

    # Alignment coordinates
    duckdb_append_int32(pw.appender, position)
    duckdb_append_int32(pw.appender, end_position)
    duckdb_append_uint8(pw.appender, mapq)
    duckdb_append_uint16(pw.appender, flag)
    duckdb_append_uint16(pw.appender, alignment_length)

    # Paired-end info
    duckdb_append_int32(pw.appender, template_length)
    duckdb_append_int32(pw.appender, mate_ref_id)
    duckdb_append_int32(pw.appender, mate_position)

    # Alignment quality tags
    if alignment_score >= 0:
        duckdb_append_int32(pw.appender, alignment_score)
    else:
        duckdb_append_null(pw.appender)

    if xs_score >= 0:
        duckdb_append_int32(pw.appender, xs_score)
    else:
        duckdb_append_null(pw.appender)

    if edit_distance >= 0:
        duckdb_append_uint16(pw.appender, <uint16_t>edit_distance)
    else:
        duckdb_append_null(pw.appender)

    if num_mismatches >= 0:
        duckdb_append_uint8(pw.appender, <uint8_t>num_mismatches)
    else:
        duckdb_append_null(pw.appender)

    if num_gap_opens >= 0:
        duckdb_append_uint8(pw.appender, <uint8_t>num_gap_opens)
    else:
        duckdb_append_null(pw.appender)

    if num_gap_extensions >= 0:
        duckdb_append_uint8(pw.appender, <uint8_t>num_gap_extensions)
    else:
        duckdb_append_null(pw.appender)

    if md_string[0] != 0:
        duckdb_append_varchar(pw.appender, md_string)
    else:
        duckdb_append_null(pw.appender)

    if ani >= 0:
        duckdb_append_float(pw.appender, ani)
    else:
        duckdb_append_null(pw.appender)

    # Ancient DNA damage
    if pmd_score >= 0:
        duckdb_append_float(pw.appender, pmd_score)
    else:
        duckdb_append_null(pw.appender)

    # bam-filter pipeline state
    if zs_score >= 0:
        duckdb_append_float(pw.appender, zs_score)
    else:
        duckdb_append_null(pw.appender)

    if zp_posterior >= 0:
        duckdb_append_float(pw.appender, zp_posterior)
    else:
        duckdb_append_null(pw.appender)

    if lca_taxid >= 0:
        duckdb_append_int32(pw.appender, lca_taxid)
    else:
        duckdb_append_null(pw.appender)

    if reassigned_ref_id >= 0:
        duckdb_append_int32(pw.appender, reassigned_ref_id)
    else:
        duckdb_append_null(pw.appender)

    # Filter status (BOOLEAN in DuckDB)
    if filter_passed:
        duckdb_append_uint8(pw.appender, 1)
    else:
        duckdb_append_uint8(pw.appender, 0)

    # Read group
    if read_group[0] != 0:
        duckdb_append_varchar(pw.appender, read_group)
    else:
        duckdb_append_null(pw.appender)

    # Sequence data
    duckdb_append_varchar(pw.appender, cigar_buf)
    duckdb_append_blob(pw.appender, seq, (alignment_length + 1) >> 1)
    duckdb_append_blob(pw.appender, qual, alignment_length)

    # Check if end_row succeeds
    if duckdb_appender_end_row(pw.appender) == DuckDBError:
        # Log error but continue - don't fail entire conversion
        bf_nogil_logf_verbose(1, "STREAM", "ERROR: Failed to append row for read_id %llu\n", read_id)
        return -1

    pw.records_written += 1

    # Flush every 10K records to ensure data is written to disk
    if pw.records_written % 10000 == 0:
        if duckdb_appender_flush(pw.appender) == DuckDBError:
            bf_nogil_logf_verbose(1, "STREAM", "ERROR: Failed to flush appender at %llu records\n", pw.records_written)
            return -1

    return 0


cdef int partition_writer_flush(PartitionWriter* pw) except -1 nogil:
    """Flush appender to disk."""
    if pw != NULL and pw.is_open:
        duckdb_appender_flush(pw.appender)
    return 0


cdef void partition_writer_close(PartitionWriter* pw) nogil:
    """Close partition writer and write to Parquet."""
    cdef char parquet_path[512]
    cdef char query[2048]

    if pw != NULL and pw.is_open:
        duckdb_appender_flush(pw.appender)
        duckdb_appender_destroy(&pw.appender)

        # Export to Parquet
        snprintf(parquet_path, sizeof(parquet_path), "%s.parquet",
                pw.filepath)
        snprintf(query, sizeof(query),
                "COPY (SELECT * FROM %s) TO '%s' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000)",
                pw.tablename, parquet_path)

        duckdb_query(pw.conn, query, NULL)

        duckdb_disconnect(&pw.conn)
        duckdb_close(&pw.db)
        pw.is_open = False


# ============================================================================
# Main Conversion with Parallel Processing
# ============================================================================

def convert_sam_bam_to_parquet_streaming(
    str input_file,
    str output_dir,
    int num_partitions=256,
    int batch_size=100000,
    int num_threads=4,
    bint write_by_reference=True,
    bint write_by_read=True,
):
    """
    Streaming Parquet conversion with two-pass processing.

    To avoid having 512+ open DuckDB connections (256 ref + 256 read), we process
    in two passes if both write_by_reference and write_by_read are requested:
    - Pass 1: Write by_reference partitions
    - Pass 2: Write by_read partitions

    This keeps the number of open databases manageable while maintaining correctness.
    """
    import time
    import sys
    import threading
    from pathlib import Path

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    if write_by_reference:
        (Path(output_dir) / "by_reference").mkdir(exist_ok=True)
    if write_by_read:
        (Path(output_dir) / "by_read").mkdir(exist_ok=True)

    # Use two-pass processing if both outputs requested to avoid too many open files
    if write_by_reference and write_by_read:
        print(f"Starting TWO-PASS streaming conversion ({num_partitions} partitions)...", flush=True)
        print(f"  Pass 1: by_reference partitions", flush=True)
        stats1 = _convert_sam_bam_to_parquet_single_pass(
            input_file, output_dir, num_partitions, batch_size, num_threads,
            write_by_reference=True, write_by_read=False)

        print(f"\n  Pass 2: by_read partitions", flush=True)
        stats2 = _convert_sam_bam_to_parquet_single_pass(
            input_file, output_dir, num_partitions, batch_size, num_threads,
            write_by_reference=False, write_by_read=True)

        # Combine stats
        return {
            'total_records': stats1['total_records'],
            'processing_time_seconds': stats1['processing_time_seconds'] + stats2['processing_time_seconds'],
        }
    else:
        print(f"Starting streaming conversion ({num_partitions} partitions, batched writes)...", flush=True)
        return _convert_sam_bam_to_parquet_single_pass(
            input_file, output_dir, num_partitions, batch_size, num_threads,
            write_by_reference, write_by_read)


cdef _convert_sam_bam_to_parquet_single_pass(
    str input_file,
    str output_dir,
    int num_partitions,
    int batch_size,
    int num_threads,
    bint write_by_reference,
    bint write_by_read,
):
    """Single-pass conversion (internal function)."""
    import time

    cdef bytes input_file_bytes = input_file.encode('utf-8')
    cdef bytes output_dir_bytes = output_dir.encode('utf-8')
    cdef const char* input_cstr = input_file_bytes
    cdef const char* output_cstr = output_dir_bytes

    cdef double start_time = time.perf_counter()

    # C-level implementation starts here - no GIL needed
    cdef:
        with nogil:
            _convert_streaming_impl(
                input_cstr,
                output_cstr,
                num_partitions,
                batch_size,
                num_threads,
                write_by_reference,
                write_by_read,
                &total_records,
                &processing_time
            )
    finally:
        running = False
        monitor_thread.join(timeout=1)

    cdef double duration = time.perf_counter() - start_time

    print(f"\nConversion complete: {total_records:,} records in {duration:.1f}s ({total_records/duration/1e6:.2f} M/sec)", flush=True)

    return {
        'total_records': total_records,
        'processing_time_seconds': processing_time,
        'duration_seconds': duration,
    }


# Batch of BAM records for parallel processing
cdef struct ReadBatch:
    bam1_t** records
    uint64_t* read_ids
    int count
    int capacity


cdef ReadBatch* read_batch_create(int capacity) nogil:
    """Create a batch buffer for records."""
    cdef ReadBatch* batch = <ReadBatch*>malloc(sizeof(ReadBatch))
    batch.records = <bam1_t**>malloc(capacity * sizeof(bam1_t*))
    batch.read_ids = <uint64_t*>malloc(capacity * sizeof(uint64_t))
    batch.count = 0
    batch.capacity = capacity

    cdef int i
    for i in range(capacity):
        batch.records[i] = bam_init1()

    return batch


cdef void read_batch_destroy(ReadBatch* batch) nogil:
    """Free batch resources."""
    cdef int i

    if batch != NULL:
        if batch.records != NULL:
            for i in range(batch.capacity):
                if batch.records[i] != NULL:
                    bam_destroy1(batch.records[i])
            free(batch.records)
        if batch.read_ids != NULL:
            free(batch.read_ids)
        free(batch)


cdef void _convert_streaming_impl(
    const char* input_file,
    const char* output_dir,
    int num_partitions,
    int batch_size,
    int num_threads,
    bint write_by_reference,
    bint write_by_read,
    uint64_t* out_total_records,
    double* out_processing_time
) nogil:
    """Streaming implementation with parallel batch processing."""
    cdef double start_time = bf_monotonic_seconds()

    # Open input
    cdef samFile* sam_fp = sam_open(input_file, "r")
    if sam_fp == NULL:
        with gil:
            raise IOError(f"Failed to open: {input_file.decode('utf-8')}")

    cdef bam_hdr_t* header = sam_hdr_read(sam_fp)
    if header == NULL:
        sam_close(sam_fp)
        with gil:
            raise IOError("Failed to read header")

    # Create partition writers (single-threaded to avoid DuckDB thread-safety issues)
    cdef PartitionWriter*** ref_writers_per_thread = NULL
    cdef PartitionWriter*** read_writers_per_thread = NULL
    cdef int t, i
    cdef char ref_dir[512]
    cdef char read_dir[512]
    cdef int read_result
    cdef uint64_t total_reads
    cdef double last_log
    cdef double current_time
    cdef uint32_t ref_partition, read_partition
    cdef ReadBatch* batch

    # Allocate for single thread only (thread 0)
    if write_by_reference:
        ref_writers_per_thread = <PartitionWriter***>calloc(1, sizeof(PartitionWriter**))
        ref_writers_per_thread[0] = <PartitionWriter**>calloc(num_partitions, sizeof(PartitionWriter*))

    if write_by_read:
        read_writers_per_thread = <PartitionWriter***>calloc(1, sizeof(PartitionWriter**))
        read_writers_per_thread[0] = <PartitionWriter**>calloc(num_partitions, sizeof(PartitionWriter*))

    snprintf(ref_dir, sizeof(ref_dir), "%s/by_reference", output_dir)
    snprintf(read_dir, sizeof(read_dir), "%s/by_read", output_dir)

    # Don't pre-create writers - create on-demand to avoid hitting system limits
    # with 256+ partitions (would be 512+ open database files)
    bf_nogil_logf_verbose(2, "STREAM", "Starting streaming conversion (single-threaded DuckDB writes, on-demand writer creation)\n")

    # Batch processing
    batch = read_batch_create(batch_size)
    total_reads = 0
    last_log = start_time

    # Read and process in batches
    while True:
        # Fill batch (sequential read - HTSlib limitation)
        batch.count = 0
        while batch.count < batch_size:
            read_result = sam_read1(sam_fp, header, batch.records[batch.count])
            if read_result < 0:
                break

            if batch.records[batch.count].core.tid >= 0:
                batch.read_ids[batch.count] = total_reads + batch.count
                batch.count += 1

        if batch.count == 0:
            break

        # Process batch serially (DuckDB appenders are not thread-safe)
        for i in range(batch.count):
            # Hash partitioning
            ref_partition = <uint32_t>(batch.records[i].core.tid % num_partitions)
            read_partition = <uint32_t>(batch.read_ids[i] % num_partitions)

            # Create writers on-demand
            if write_by_reference:
                if ref_writers_per_thread[0][ref_partition] == NULL:
                    ref_writers_per_thread[0][ref_partition] = partition_writer_create(
                        ref_dir, ref_partition, "by_ref", True, 0)
                    if ref_writers_per_thread[0][ref_partition] == NULL:
                        bf_nogil_logf_verbose(1, "STREAM", "ERROR: Failed to create ref writer for partition %u\n", ref_partition)
                        continue

                if partition_writer_add_record(ref_writers_per_thread[0][ref_partition],
                                           batch.records[i], header, batch.read_ids[i]) < 0:
                    bf_nogil_logf_verbose(1, "STREAM", "Failed to write record %llu to ref partition %u\n",
                                         batch.read_ids[i], ref_partition)

            if write_by_read:
                if read_writers_per_thread[0][read_partition] == NULL:
                    read_writers_per_thread[0][read_partition] = partition_writer_create(
                        read_dir, read_partition, "by_read", False, 0)
                    if read_writers_per_thread[0][read_partition] == NULL:
                        bf_nogil_logf_verbose(1, "STREAM", "ERROR: Failed to create read writer for partition %u\n", read_partition)
                        continue

                if partition_writer_add_record(read_writers_per_thread[0][read_partition],
                                           batch.records[i], header, batch.read_ids[i]) < 0:
                    bf_nogil_logf_verbose(1, "STREAM", "Failed to write record %llu to read partition %u\n",
                                         batch.read_ids[i], read_partition)

        total_reads += batch.count

        # Progress logging
        current_time = bf_monotonic_seconds()
        if current_time - last_log >= 10.0:
            bf_nogil_logf_verbose(2, "STREAM",
                "Processed %llu records (%.1f M/sec)\n",
                total_reads, total_reads / (current_time - start_time) / 1000000.0)
            last_log = current_time

    # Close all writers and export to Parquet (only thread 0 is used now)
    if write_by_reference:
        for i in range(num_partitions):
            if ref_writers_per_thread[0][i] != NULL:
                partition_writer_close(ref_writers_per_thread[0][i])
                free(ref_writers_per_thread[0][i])
        free(ref_writers_per_thread[0])
        free(ref_writers_per_thread)

    if write_by_read:
        for i in range(num_partitions):
            if read_writers_per_thread[0][i] != NULL:
                partition_writer_close(read_writers_per_thread[0][i])
                free(read_writers_per_thread[0][i])
        free(read_writers_per_thread[0])
        free(read_writers_per_thread)

    # Cleanup
    read_batch_destroy(batch)
    bam_hdr_destroy(header)
    sam_close(sam_fp)

    cdef double end_time = bf_monotonic_seconds()
    out_total_records[0] = total_reads
    out_processing_time[0] = end_time - start_time

    bf_nogil_logf_verbose(2, "STREAM",
        "Complete: %llu records in %.1f sec (%.1f M/sec)\n",
        total_reads, end_time - start_time,
        total_reads / (end_time - start_time) / 1000000.0)
