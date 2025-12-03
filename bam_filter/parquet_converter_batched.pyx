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
Batched Parquet converter with smart buffering.

Key improvements:
- Accumulates 100K records per partition in memory (columnar format)
- Bulk appends to DuckDB (minimize API calls)
- Single-pass processing for both by_reference and by_read
- Lazy partition file creation
- Target: 10-100x faster than row-by-row appends
"""

from libc.stdlib cimport malloc, free, calloc, realloc
from libc.string cimport memcpy, strlen, strcpy, strcmp, memset, strdup
from libc.stdint cimport uint8_t, uint16_t, uint32_t, int32_t, uint64_t, int8_t, int16_t
from libc.stdio cimport snprintf, FILE, fopen, fclose, fprintf
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


# ============================================================================
# Tag Extraction (same as before)
# ============================================================================

cdef inline int32_t extract_tag_int(bam1_t* record, const char* tag_name) nogil:
    """Extract integer tag (AS, NM, etc.)."""
    cdef uint8_t* aux = bam_get_aux(record)
    cdef uint8_t* end = record.data + record.l_data
    cdef uint8_t tag_type

    while aux < end - 2:
        if aux[0] == tag_name[0] and aux[1] == tag_name[1]:
            aux += 2
            if aux[0] == b'C':
                return <int32_t>(<uint8_t*>(aux + 1))[0]
            elif aux[0] == b'c':
                return <int32_t>(<int8_t*>(aux + 1))[0]
            elif aux[0] == b'S':
                return <int32_t>(<uint16_t*>(aux + 1))[0]
            elif aux[0] == b's':
                return <int32_t>(<int16_t*>(aux + 1))[0]
            elif aux[0] == b'I':
                return <int32_t>(<uint32_t*>(aux + 1))[0]
            elif aux[0] == b'i':
                return (<int32_t*>(aux + 1))[0]
            return -1

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


cdef inline float extract_tag_float(bam1_t* record, const char* tag_name) nogil:
    """Extract float tag (PMD, ZS, ZP, etc.)."""
    cdef uint8_t* aux = bam_get_aux(record)
    cdef uint8_t* end = record.data + record.l_data
    cdef uint8_t tag_type

    if aux == NULL or aux >= end - 2:
        return -1.0

    while aux < end - 2:
        if aux[0] == tag_name[0] and aux[1] == tag_name[1]:
            aux += 2
            tag_type = aux[0]
            aux += 1

            if tag_type == b'f':
                return (<float*>aux)[0]
            elif tag_type == b'd':
                return <float>((<double*>aux)[0])
            elif tag_type == b'i':
                return <float>((<int32_t*>aux)[0])
            elif tag_type == b'I':
                return <float>((<uint32_t*>aux)[0])
            return -1.0

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
            aux += 1
        else:
            break

    return -1.0


cdef inline int extract_tag_string(bam1_t* record, const char* tag_name, char* buffer, int max_len) nogil:
    """Extract string tag (MD, RG, etc.)."""
    cdef uint8_t* aux = bam_get_aux(record)
    cdef uint8_t* end = record.data + record.l_data
    cdef int len_written = 0
    cdef uint8_t tag_type

    buffer[0] = 0

    while aux < end - 2:
        if aux[0] == tag_name[0] and aux[1] == tag_name[1]:
            aux += 2
            if aux[0] == b'Z':
                aux += 1
                while aux < end and aux[0] != 0 and len_written < max_len - 1:
                    buffer[len_written] = <char>aux[0]
                    len_written += 1
                    aux += 1
                buffer[len_written] = 0
                return len_written
            return 0

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


# ============================================================================
# Columnar Batch Buffer - stores records in columnar format
# ============================================================================

cdef struct PartitionBatch:
    # Capacity and current count
    uint32_t capacity
    uint32_t count
    uint32_t partition_id

    # Columnar arrays (parallel arrays)
    uint64_t* read_ids
    char** read_names
    uint32_t* ref_ids
    char** ref_names
    int32_t* positions
    int32_t* end_positions
    uint8_t* mapqs
    uint16_t* flags
    uint16_t* alignment_lengths
    int32_t* template_lengths
    int32_t* mate_ref_ids
    int32_t* mate_positions

    # Tags
    int32_t* alignment_scores
    int32_t* xs_scores
    int32_t* edit_distances
    int32_t* num_mismatches
    int32_t* num_gap_opens
    int32_t* num_gap_extensions
    char** md_strings
    float* anis
    float* pmd_scores
    float* zs_scores
    float* zp_posteriors
    int32_t* lca_taxids
    int32_t* reassigned_ref_ids
    uint8_t* filter_passed
    char** read_groups
    char** cigars
    uint8_t** sequences
    uint16_t* seq_lengths  # For blob data
    uint8_t** qualities


cdef PartitionBatch* partition_batch_create(uint32_t partition_id, uint32_t capacity) nogil:
    """Create a columnar batch buffer."""
    cdef PartitionBatch* batch = <PartitionBatch*>calloc(1, sizeof(PartitionBatch))
    if batch == NULL:
        return NULL

    batch.partition_id = partition_id
    batch.capacity = capacity
    batch.count = 0

    # Allocate columnar arrays
    batch.read_ids = <uint64_t*>malloc(capacity * sizeof(uint64_t))
    batch.read_names = <char**>malloc(capacity * sizeof(char*))
    batch.ref_ids = <uint32_t*>malloc(capacity * sizeof(uint32_t))
    batch.ref_names = <char**>malloc(capacity * sizeof(char*))
    batch.positions = <int32_t*>malloc(capacity * sizeof(int32_t))
    batch.end_positions = <int32_t*>malloc(capacity * sizeof(int32_t))
    batch.mapqs = <uint8_t*>malloc(capacity * sizeof(uint8_t))
    batch.flags = <uint16_t*>malloc(capacity * sizeof(uint16_t))
    batch.alignment_lengths = <uint16_t*>malloc(capacity * sizeof(uint16_t))
    batch.template_lengths = <int32_t*>malloc(capacity * sizeof(int32_t))
    batch.mate_ref_ids = <int32_t*>malloc(capacity * sizeof(int32_t))
    batch.mate_positions = <int32_t*>malloc(capacity * sizeof(int32_t))

    batch.alignment_scores = <int32_t*>malloc(capacity * sizeof(int32_t))
    batch.xs_scores = <int32_t*>malloc(capacity * sizeof(int32_t))
    batch.edit_distances = <int32_t*>malloc(capacity * sizeof(int32_t))
    batch.num_mismatches = <int32_t*>malloc(capacity * sizeof(int32_t))
    batch.num_gap_opens = <int32_t*>malloc(capacity * sizeof(int32_t))
    batch.num_gap_extensions = <int32_t*>malloc(capacity * sizeof(int32_t))
    batch.md_strings = <char**>malloc(capacity * sizeof(char*))
    batch.anis = <float*>malloc(capacity * sizeof(float))
    batch.pmd_scores = <float*>malloc(capacity * sizeof(float))
    batch.zs_scores = <float*>malloc(capacity * sizeof(float))
    batch.zp_posteriors = <float*>malloc(capacity * sizeof(float))
    batch.lca_taxids = <int32_t*>malloc(capacity * sizeof(int32_t))
    batch.reassigned_ref_ids = <int32_t*>malloc(capacity * sizeof(int32_t))
    batch.filter_passed = <uint8_t*>malloc(capacity * sizeof(uint8_t))
    batch.read_groups = <char**>malloc(capacity * sizeof(char*))
    batch.cigars = <char**>malloc(capacity * sizeof(char*))
    batch.sequences = <uint8_t**>malloc(capacity * sizeof(uint8_t*))
    batch.seq_lengths = <uint16_t*>malloc(capacity * sizeof(uint16_t))
    batch.qualities = <uint8_t**>malloc(capacity * sizeof(uint8_t*))

    return batch


cdef void partition_batch_destroy(PartitionBatch* batch) nogil:
    """Free batch resources."""
    cdef uint32_t i

    if batch == NULL:
        return

    # Free string arrays
    if batch.read_names != NULL:
        for i in range(batch.count):
            if batch.read_names[i] != NULL:
                free(batch.read_names[i])
        free(batch.read_names)

    if batch.ref_names != NULL:
        for i in range(batch.count):
            if batch.ref_names[i] != NULL:
                free(batch.ref_names[i])
        free(batch.ref_names)

    if batch.md_strings != NULL:
        for i in range(batch.count):
            if batch.md_strings[i] != NULL:
                free(batch.md_strings[i])
        free(batch.md_strings)

    if batch.read_groups != NULL:
        for i in range(batch.count):
            if batch.read_groups[i] != NULL:
                free(batch.read_groups[i])
        free(batch.read_groups)

    if batch.cigars != NULL:
        for i in range(batch.count):
            if batch.cigars[i] != NULL:
                free(batch.cigars[i])
        free(batch.cigars)

    if batch.sequences != NULL:
        for i in range(batch.count):
            if batch.sequences[i] != NULL:
                free(batch.sequences[i])
        free(batch.sequences)

    if batch.qualities != NULL:
        for i in range(batch.count):
            if batch.qualities[i] != NULL:
                free(batch.qualities[i])
        free(batch.qualities)

    # Free scalar arrays
    free(batch.read_ids)
    free(batch.ref_ids)
    free(batch.positions)
    free(batch.end_positions)
    free(batch.mapqs)
    free(batch.flags)
    free(batch.alignment_lengths)
    free(batch.template_lengths)
    free(batch.mate_ref_ids)
    free(batch.mate_positions)
    free(batch.alignment_scores)
    free(batch.xs_scores)
    free(batch.edit_distances)
    free(batch.num_mismatches)
    free(batch.num_gap_opens)
    free(batch.num_gap_extensions)
    free(batch.anis)
    free(batch.pmd_scores)
    free(batch.zs_scores)
    free(batch.zp_posteriors)
    free(batch.lca_taxids)
    free(batch.reassigned_ref_ids)
    free(batch.filter_passed)
    free(batch.seq_lengths)

    free(batch)


cdef int partition_batch_add_record(PartitionBatch* batch, bam1_t* record, bam_hdr_t* header,
                                     uint64_t read_id) except -1 nogil:
    """Add record to batch (columnar storage)."""
    if batch.count >= batch.capacity:
        return -1  # Batch full

    cdef uint32_t idx = batch.count

    # Extract basic fields
    batch.read_ids[idx] = read_id
    batch.ref_ids[idx] = <uint32_t>record.core.tid
    batch.positions[idx] = record.core.pos
    batch.end_positions[idx] = bam_endpos(record)
    batch.mapqs[idx] = record.core.qual
    batch.flags[idx] = record.core.flag
    batch.alignment_lengths[idx] = <uint16_t>record.core.l_qseq
    batch.template_lengths[idx] = record.core.isize
    batch.mate_ref_ids[idx] = record.core.mtid
    batch.mate_positions[idx] = record.core.mpos

    # Copy read name
    cdef char* read_name = bam_get_qname(record)
    batch.read_names[idx] = strdup(read_name)

    # Copy reference name
    cdef const char* ref_name
    if batch.ref_ids[idx] >= 0 and batch.ref_ids[idx] < header.n_targets and header.target_name != NULL:
        ref_name = header.target_name[batch.ref_ids[idx]]
        if ref_name != NULL:
            batch.ref_names[idx] = strdup(ref_name)
        else:
            batch.ref_names[idx] = strdup("*")
    else:
        batch.ref_names[idx] = strdup("*")

    # Extract tags
    batch.alignment_scores[idx] = extract_tag_int(record, "AS")
    batch.xs_scores[idx] = extract_tag_int(record, "XS")
    batch.edit_distances[idx] = extract_tag_int(record, "NM")
    batch.num_mismatches[idx] = extract_tag_int(record, "XM")
    batch.num_gap_opens[idx] = extract_tag_int(record, "XO")
    batch.num_gap_extensions[idx] = extract_tag_int(record, "XG")
    batch.pmd_scores[idx] = extract_tag_float(record, "PM")
    if batch.pmd_scores[idx] < 0:
        batch.pmd_scores[idx] = extract_tag_float(record, "PMD")
    batch.zs_scores[idx] = extract_tag_float(record, "ZS")
    batch.zp_posteriors[idx] = extract_tag_float(record, "ZP")
    batch.lca_taxids[idx] = extract_tag_int(record, "ZT")
    batch.reassigned_ref_ids[idx] = extract_tag_int(record, "ZR")

    cdef int32_t filter_flag = extract_tag_int(record, "ZF")
    batch.filter_passed[idx] = 1 if (filter_flag != 0 or filter_flag < 0) else 0

    # Calculate ANI
    if batch.edit_distances[idx] >= 0 and batch.alignment_lengths[idx] > 0:
        batch.anis[idx] = (1.0 - (<float>batch.edit_distances[idx] / <float>batch.alignment_lengths[idx])) * 100.0
    else:
        batch.anis[idx] = -1.0

    # Extract MD and RG tags
    cdef char md_buf[1024]
    cdef char rg_buf[256]
    extract_tag_string(record, "MD", md_buf, sizeof(md_buf))
    extract_tag_string(record, "RG", rg_buf, sizeof(rg_buf))
    batch.md_strings[idx] = strdup(md_buf)
    batch.read_groups[idx] = strdup(rg_buf)

    # Build CIGAR string
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
        if op < 9:
            op_char = cigar_ops[op]
        else:
            op_char = '?'
        pos += snprintf(cigar_buf + pos, sizeof(cigar_buf) - pos, "%u%c",
                       op_len, op_char)
        if pos >= sizeof(cigar_buf) - 1:
            break
    cigar_buf[pos] = 0
    batch.cigars[idx] = strdup(cigar_buf)

    # Copy sequence and quality
    cdef uint8_t* seq = bam_get_seq(record)
    cdef uint8_t* qual = bam_get_qual(record)
    cdef uint16_t seq_len = batch.alignment_lengths[idx]
    batch.seq_lengths[idx] = seq_len

    batch.sequences[idx] = <uint8_t*>malloc((seq_len + 1) >> 1)
    memcpy(batch.sequences[idx], seq, (seq_len + 1) >> 1)

    batch.qualities[idx] = <uint8_t*>malloc(seq_len)
    memcpy(batch.qualities[idx], qual, seq_len)

    batch.count += 1
    return 0


# ============================================================================
# Partition Writer - opens DuckDB file and writes batches
# ============================================================================

cdef struct PartitionWriter:
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
                                               bint is_reference_table) except NULL nogil:
    """Create partition writer."""
    cdef PartitionWriter* pw = <PartitionWriter*>calloc(1, sizeof(PartitionWriter))
    if pw == NULL:
        return NULL

    pw.partition_id = partition_id
    pw.records_written = 0
    pw.is_open = False

    # Create file path
    if is_reference_table:
        snprintf(pw.filepath, sizeof(pw.filepath),
                "%s/by_reference/by_reference_p%04d.duckdb", output_dir, partition_id)
    else:
        snprintf(pw.filepath, sizeof(pw.filepath),
                "%s/by_read/by_read_p%04d.duckdb", output_dir, partition_id)

    snprintf(pw.tablename, sizeof(pw.tablename), "alignments")

    # Open DuckDB
    if duckdb_open(pw.filepath, &pw.db) == DuckDBError:
        free(pw)
        return NULL

    if duckdb_connect(pw.db, &pw.conn) == DuckDBError:
        duckdb_close(&pw.db)
        free(pw)
        return NULL

    # Create table
    cdef char query[4096]
    snprintf(query, sizeof(query),
            "CREATE TABLE %s ("
            "read_id UBIGINT, "
            "read_name VARCHAR, "
            "ref_id UINTEGER, "
            "ref_name VARCHAR, "
            "position INTEGER, "
            "end_position INTEGER, "
            "mapq UTINYINT, "
            "flag USMALLINT, "
            "alignment_length USMALLINT, "
            "template_length INTEGER, "
            "mate_ref_id INTEGER, "
            "mate_position INTEGER, "
            "alignment_score INTEGER, "
            "xs_score INTEGER, "
            "edit_distance USMALLINT, "
            "num_mismatches UTINYINT, "
            "num_gap_opens UTINYINT, "
            "num_gap_extensions UTINYINT, "
            "md_string VARCHAR, "
            "ani FLOAT, "
            "pmd_score FLOAT, "
            "zs_score FLOAT, "
            "zp_posterior FLOAT, "
            "lca_taxid INTEGER, "
            "reassigned_ref_id INTEGER, "
            "filter_passed BOOLEAN, "
            "read_group VARCHAR, "
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


cdef int partition_writer_write_batch(PartitionWriter* pw, PartitionBatch* batch) except -1 nogil:
    """Write entire batch to DuckDB using appender."""
    if pw == NULL or not pw.is_open or batch == NULL:
        return -1

    cdef uint32_t i

    # Bulk append all rows
    for i in range(batch.count):
        duckdb_appender_begin_row(pw.appender)

        # Core fields
        duckdb_append_uint64(pw.appender, batch.read_ids[i])
        duckdb_append_varchar(pw.appender, batch.read_names[i])
        duckdb_append_uint32(pw.appender, batch.ref_ids[i])
        duckdb_append_varchar(pw.appender, batch.ref_names[i])
        duckdb_append_int32(pw.appender, batch.positions[i])
        duckdb_append_int32(pw.appender, batch.end_positions[i])
        duckdb_append_uint8(pw.appender, batch.mapqs[i])
        duckdb_append_uint16(pw.appender, batch.flags[i])
        duckdb_append_uint16(pw.appender, batch.alignment_lengths[i])
        duckdb_append_int32(pw.appender, batch.template_lengths[i])
        duckdb_append_int32(pw.appender, batch.mate_ref_ids[i])
        duckdb_append_int32(pw.appender, batch.mate_positions[i])

        # Tags with NULL handling
        if batch.alignment_scores[i] >= 0:
            duckdb_append_int32(pw.appender, batch.alignment_scores[i])
        else:
            duckdb_append_null(pw.appender)

        if batch.xs_scores[i] >= 0:
            duckdb_append_int32(pw.appender, batch.xs_scores[i])
        else:
            duckdb_append_null(pw.appender)

        if batch.edit_distances[i] >= 0:
            duckdb_append_uint16(pw.appender, <uint16_t>batch.edit_distances[i])
        else:
            duckdb_append_null(pw.appender)

        if batch.num_mismatches[i] >= 0:
            duckdb_append_uint8(pw.appender, <uint8_t>batch.num_mismatches[i])
        else:
            duckdb_append_null(pw.appender)

        if batch.num_gap_opens[i] >= 0:
            duckdb_append_uint8(pw.appender, <uint8_t>batch.num_gap_opens[i])
        else:
            duckdb_append_null(pw.appender)

        if batch.num_gap_extensions[i] >= 0:
            duckdb_append_uint8(pw.appender, <uint8_t>batch.num_gap_extensions[i])
        else:
            duckdb_append_null(pw.appender)

        if batch.md_strings[i][0] != 0:
            duckdb_append_varchar(pw.appender, batch.md_strings[i])
        else:
            duckdb_append_null(pw.appender)

        if batch.anis[i] >= 0:
            duckdb_append_float(pw.appender, batch.anis[i])
        else:
            duckdb_append_null(pw.appender)

        if batch.pmd_scores[i] >= 0:
            duckdb_append_float(pw.appender, batch.pmd_scores[i])
        else:
            duckdb_append_null(pw.appender)

        if batch.zs_scores[i] >= 0:
            duckdb_append_float(pw.appender, batch.zs_scores[i])
        else:
            duckdb_append_null(pw.appender)

        if batch.zp_posteriors[i] >= 0:
            duckdb_append_float(pw.appender, batch.zp_posteriors[i])
        else:
            duckdb_append_null(pw.appender)

        if batch.lca_taxids[i] >= 0:
            duckdb_append_int32(pw.appender, batch.lca_taxids[i])
        else:
            duckdb_append_null(pw.appender)

        if batch.reassigned_ref_ids[i] >= 0:
            duckdb_append_int32(pw.appender, batch.reassigned_ref_ids[i])
        else:
            duckdb_append_null(pw.appender)

        duckdb_append_uint8(pw.appender, batch.filter_passed[i])

        if batch.read_groups[i][0] != 0:
            duckdb_append_varchar(pw.appender, batch.read_groups[i])
        else:
            duckdb_append_null(pw.appender)

        # Sequence data
        duckdb_append_varchar(pw.appender, batch.cigars[i])
        duckdb_append_blob(pw.appender, batch.sequences[i], (batch.seq_lengths[i] + 1) >> 1)
        duckdb_append_blob(pw.appender, batch.qualities[i], batch.seq_lengths[i])

        if duckdb_appender_end_row(pw.appender) == DuckDBError:
            bf_nogil_logf_verbose(1, "BATCH", "ERROR: Failed to append row %u\n", i)
            return -1

    # Flush after batch
    if duckdb_appender_flush(pw.appender) == DuckDBError:
        bf_nogil_logf_verbose(1, "BATCH", "ERROR: Failed to flush batch\n")
        return -1

    pw.records_written += batch.count
    return 0


cdef void partition_writer_close(PartitionWriter* pw) nogil:
    """Close writer and export to Parquet."""
    cdef char parquet_path[512]
    cdef char query[2048]

    if pw != NULL and pw.is_open:
        duckdb_appender_flush(pw.appender)
        duckdb_appender_destroy(&pw.appender)

        # Export to Parquet
        snprintf(parquet_path, sizeof(parquet_path), "%s.parquet", pw.filepath)
        snprintf(query, sizeof(query),
                "COPY (SELECT * FROM %s) TO '%s' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000)",
                pw.tablename, parquet_path)

        duckdb_query(pw.conn, query, NULL)

        duckdb_disconnect(&pw.conn)
        duckdb_close(&pw.db)
        pw.is_open = False


# ============================================================================
# Main Conversion Function
# ============================================================================

def convert_sam_bam_to_parquet_batched(
    str input_file,
    str output_dir,
    int num_partitions=128,
    int batch_size=100000,
    bint write_by_reference=True,
    bint write_by_read=True,
):
    """
    Batched Parquet conversion with smart buffering.

    Accumulates records in memory per partition, then bulk writes when buffers full.
    Single-pass processing for both by_reference and by_read outputs.
    """
    import time
    from pathlib import Path

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    if write_by_reference:
        (Path(output_dir) / "by_reference").mkdir(exist_ok=True)
    if write_by_read:
        (Path(output_dir) / "by_read").mkdir(exist_ok=True)

    cdef bytes input_file_bytes = input_file.encode('utf-8')
    cdef bytes output_dir_bytes = output_dir.encode('utf-8')
    cdef const char* input_cstr = input_file_bytes
    cdef const char* output_cstr = output_dir_bytes

    cdef double start_time = time.perf_counter()
    cdef uint64_t total_records = 0

    with nogil:
        total_records = _convert_batched_impl(
            input_cstr,
            output_cstr,
            num_partitions,
            batch_size,
            write_by_reference,
            write_by_read
        )

    cdef double duration = time.perf_counter() - start_time

    print(f"\nConversion complete: {total_records:,} records in {duration:.1f}s ({total_records/duration/1e6:.2f} M/sec)")

    return {
        'total_records': total_records,
        'processing_time_seconds': duration,
    }


cdef uint64_t _convert_batched_impl(
    const char* input_file,
    const char* output_dir,
    int num_partitions,
    int batch_size,
    bint write_by_reference,
    bint write_by_read,
) nogil:
    """Batched conversion implementation."""
    cdef double start_time = bf_monotonic_seconds()

    # Open input
    cdef samFile* sam_fp = sam_open(input_file, "r")
    if sam_fp == NULL:
        with gil:
            raise IOError(f"Failed to open input file")

    cdef bam_hdr_t* header = sam_hdr_read(sam_fp)
    if header == NULL:
        sam_close(sam_fp)
        with gil:
            raise IOError("Failed to read header")

    # Create partition batches (lazy allocation)
    cdef PartitionBatch** ref_batches = NULL
    cdef PartitionBatch** read_batches = NULL

    if write_by_reference:
        ref_batches = <PartitionBatch**>calloc(num_partitions, sizeof(PartitionBatch*))

    if write_by_read:
        read_batches = <PartitionBatch**>calloc(num_partitions, sizeof(PartitionBatch*))

    # Create partition writers (lazy allocation)
    cdef PartitionWriter** ref_writers = NULL
    cdef PartitionWriter** read_writers = NULL

    if write_by_reference:
        ref_writers = <PartitionWriter**>calloc(num_partitions, sizeof(PartitionWriter*))

    if write_by_read:
        read_writers = <PartitionWriter**>calloc(num_partitions, sizeof(PartitionWriter*))

    # Process SAM file
    cdef bam1_t* record = bam_init1()
    cdef uint64_t total_reads = 0
    cdef int read_result
    cdef uint32_t ref_partition, read_partition
    cdef double last_log = start_time
    cdef double current_time

    bf_nogil_logf_verbose(2, "BATCH", "Starting batched conversion\n")

    while True:
        read_result = sam_read1(sam_fp, header, record)
        if read_result < 0:
            break

        if record.core.tid < 0:
            continue

        # Determine partitions
        ref_partition = <uint32_t>(record.core.tid % num_partitions)
        read_partition = <uint32_t>(total_reads % num_partitions)

        # Add to by_reference batch
        if write_by_reference:
            # Create batch if needed
            if ref_batches[ref_partition] == NULL:
                ref_batches[ref_partition] = partition_batch_create(ref_partition, batch_size)

            # Add record
            if partition_batch_add_record(ref_batches[ref_partition], record, header, total_reads) < 0:
                # Batch full - flush it
                if ref_writers[ref_partition] == NULL:
                    ref_writers[ref_partition] = partition_writer_create(output_dir, ref_partition, True)

                partition_writer_write_batch(ref_writers[ref_partition], ref_batches[ref_partition])

                # Clear batch
                partition_batch_destroy(ref_batches[ref_partition])
                ref_batches[ref_partition] = partition_batch_create(ref_partition, batch_size)
                partition_batch_add_record(ref_batches[ref_partition], record, header, total_reads)

        # Add to by_read batch
        if write_by_read:
            if read_batches[read_partition] == NULL:
                read_batches[read_partition] = partition_batch_create(read_partition, batch_size)

            if partition_batch_add_record(read_batches[read_partition], record, header, total_reads) < 0:
                if read_writers[read_partition] == NULL:
                    read_writers[read_partition] = partition_writer_create(output_dir, read_partition, False)

                partition_writer_write_batch(read_writers[read_partition], read_batches[read_partition])

                partition_batch_destroy(read_batches[read_partition])
                read_batches[read_partition] = partition_batch_create(read_partition, batch_size)
                partition_batch_add_record(read_batches[read_partition], record, header, total_reads)

        total_reads += 1

        # Progress logging
        current_time = bf_monotonic_seconds()
        if current_time - last_log >= 10.0:
            bf_nogil_logf_verbose(2, "BATCH",
                "Processed %llu records (%.1f M/sec)\n",
                total_reads, total_reads / (current_time - start_time) / 1000000.0)
            last_log = current_time

    # Flush remaining batches
    cdef int i

    if write_by_reference:
        for i in range(num_partitions):
            if ref_batches[i] != NULL and ref_batches[i].count > 0:
                if ref_writers[i] == NULL:
                    ref_writers[i] = partition_writer_create(output_dir, i, True)
                partition_writer_write_batch(ref_writers[i], ref_batches[i])
                partition_batch_destroy(ref_batches[i])

            if ref_writers[i] != NULL:
                partition_writer_close(ref_writers[i])
                free(ref_writers[i])
        free(ref_batches)
        free(ref_writers)

    if write_by_read:
        for i in range(num_partitions):
            if read_batches[i] != NULL and read_batches[i].count > 0:
                if read_writers[i] == NULL:
                    read_writers[i] = partition_writer_create(output_dir, i, False)
                partition_writer_write_batch(read_writers[i], read_batches[i])
                partition_batch_destroy(read_batches[i])

            if read_writers[i] != NULL:
                partition_writer_close(read_writers[i])
                free(read_writers[i])
        free(read_batches)
        free(read_writers)

    # Cleanup
    bam_destroy1(record)
    bam_hdr_destroy(header)
    sam_close(sam_fp)

    cdef double end_time = bf_monotonic_seconds()
    bf_nogil_logf_verbose(2, "BATCH",
        "Complete: %llu records in %.1f sec (%.1f M/sec)\n",
        total_reads, end_time - start_time,
        total_reads / (end_time - start_time) / 1000000.0)

    return total_reads
