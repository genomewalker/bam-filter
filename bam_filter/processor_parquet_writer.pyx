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
# -*- coding: utf-8 -*-

"""BAM to Parquet conversion with metagenomic-scale optimization.

Handles efficient conversion of coordinate-sorted BAM files to partitioned
Parquet format optimized for DuckDB queries. Uses streaming batch processing
to maintain constant memory usage even for billion-alignment datasets.
"""

from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t, uint8_t, uint16_t
from libc.stdlib cimport malloc, free, realloc, calloc
from libc.string cimport memcpy, memset, strlen
from libc.math cimport NAN
from libc.stddef cimport size_t
from cpython.unicode cimport PyUnicode_DecodeUTF8
from cpython.bytes cimport PyBytes_FromStringAndSize
from cpython.long cimport PyLong_FromVoidPtr
from cpython.list cimport PyList_New
from cpython.ref cimport Py_DECREF, Py_INCREF
from cpython.object cimport PyObject
from cython.parallel cimport prange, parallel, threadid

# Import HTSlib types from shared type definitions
from bam_filter.processor_types cimport (
    samFile, bam1_t, bam1_core_t, sam_hdr_t, hts_idx_t, hts_itr_t,
    bam_init1, bam_destroy1, bam_endpos,
    bam_get_qname, bam_get_seq, bam_get_qual, bam_get_cigar,
    bam_aux_get, bam_aux2i, bam_seqi_wrapper,
    hts_open, hts_close, sam_hdr_read, sam_hdr_destroy,
    sam_hdr_tid2name, sam_hdr_tid2len, sam_hdr_nref, sam_hdr_str, sam_read1,
    sam_index_load, hts_idx_destroy, hts_idx_get_stat,
    sam_itr_queryi, sam_itr_next, hts_itr_destroy
)
from bam_filter.processor cimport AlignmentScoringConfig
from bam_filter.processor_md_quality cimport calculate_md_quality_score, alignment_passes_quality_filters_with_ani
from bam_filter.batch_utils cimport create_balanced_batches_greedy

cdef extern from "bam_filter/c_logging.h":
    double bf_monotonic_seconds() nogil
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil
    int bf_should_log(int level) nogil

# External declarations
cdef extern from "stdio.h" nogil:
    int snprintf(char* s, size_t n, const char* format, ...)

cdef const char* C_LOG_TAG = b"PARQUET"


cdef inline uint8_t _encode_base_char(char base) nogil:
    cdef char b = base
    if b >= c'a' and b <= c'z':
        b -= 32
    if b == c'A':
        return 1
    elif b == c'C':
        return 2
    elif b == c'G':
        return 4
    elif b == c'T':
        return 8
    else:
        return 15


cdef inline bint pack_sequence_bases(const char* seq,
                                     int32_t length,
                                     char** out_data,
                                     uint32_t* out_size) nogil:
    if length <= 0:
        out_data[0] = NULL
        out_size[0] = 0
        return True
    cdef uint32_t packed_len = <uint32_t>((length + 1) >> 1)
    cdef char* packed = <char*>malloc(packed_len)
    if not packed:
        return False
    cdef int32_t i
    cdef uint8_t code
    cdef uint32_t idx
    for i in range(length):
        code = _encode_base_char(seq[i]) & 0x0F
        idx = <uint32_t>(i >> 1)
        if (i & 1) == 0:
            packed[idx] = <char>((code << 4) & 0xF0)
        else:
            packed[idx] = <char>(packed[idx] | code)
    out_data[0] = packed
    out_size[0] = packed_len
    return True


cdef inline bint _ensure_capacity(void** data_ptr,
                                  size_t* capacity_ptr,
                                  size_t needed,
                                  size_t element_size) nogil:
    cdef size_t current = capacity_ptr[0]
    if needed <= current:
        return True
    cdef size_t new_capacity
    if current > 0:
        new_capacity = current * 2
    else:
        new_capacity = 65536
    if new_capacity < needed:
        new_capacity = needed
    cdef void* new_data = realloc(data_ptr[0], new_capacity * element_size)
    if not new_data:
        return False
    data_ptr[0] = new_data
    capacity_ptr[0] = new_capacity
    return True


cdef struct UInt32Column:
    uint32_t* data
    size_t size
    size_t capacity


cdef struct Int32Column:
    int32_t* data
    size_t size
    size_t capacity


cdef struct UInt16Column:
    uint16_t* data
    size_t size
    size_t capacity


cdef struct UInt8Column:
    uint8_t* data
    size_t size
    size_t capacity


cdef struct FloatColumn:
    float* data
    size_t size
    size_t capacity


cdef struct VarLenColumn:
    char* data
    size_t data_size
    size_t data_capacity
    int32_t* offsets
    size_t offsets_size
    size_t offsets_capacity


cdef inline void _uint32_column_init(UInt32Column* col) nogil:
    col.data = NULL
    col.size = 0
    col.capacity = 0


cdef inline void _int32_column_init(Int32Column* col) nogil:
    col.data = NULL
    col.size = 0
    col.capacity = 0


cdef inline void _uint16_column_init(UInt16Column* col) nogil:
    col.data = NULL
    col.size = 0
    col.capacity = 0


cdef inline void _uint8_column_init(UInt8Column* col) nogil:
    col.data = NULL
    col.size = 0
    col.capacity = 0


cdef inline void _float_column_init(FloatColumn* col) nogil:
    col.data = NULL
    col.size = 0
    col.capacity = 0


cdef inline void _varlen_column_init(VarLenColumn* col) nogil:
    col.data = NULL
    col.data_size = 0
    col.data_capacity = 0
    col.offsets = NULL
    col.offsets_size = 0
    col.offsets_capacity = 0


cdef inline void _uint32_column_free(UInt32Column* col) nogil:
    if col.data:
        free(col.data)
    col.data = NULL
    col.size = 0
    col.capacity = 0


cdef inline void _int32_column_free(Int32Column* col) nogil:
    if col.data:
        free(col.data)
    col.data = NULL
    col.size = 0
    col.capacity = 0


cdef inline void _uint16_column_free(UInt16Column* col) nogil:
    if col.data:
        free(col.data)
    col.data = NULL
    col.size = 0
    col.capacity = 0


cdef inline void _uint8_column_free(UInt8Column* col) nogil:
    if col.data:
        free(col.data)
    col.data = NULL
    col.size = 0
    col.capacity = 0


cdef inline void _float_column_free(FloatColumn* col) nogil:
    if col.data:
        free(col.data)
    col.data = NULL
    col.size = 0
    col.capacity = 0


cdef inline void _varlen_column_free(VarLenColumn* col) nogil:
    if col.data:
        free(col.data)
    if col.offsets:
        free(col.offsets)
    col.data = NULL
    col.data_size = 0
    col.data_capacity = 0
    col.offsets = NULL
    col.offsets_size = 0
    col.offsets_capacity = 0


cdef inline bint _uint32_column_append(UInt32Column* col, uint32_t value) nogil:
    if not _ensure_capacity(<void**>&col.data, &col.capacity, col.size + 1, sizeof(uint32_t)):
        return False
    col.data[col.size] = value
    col.size += 1
    return True


cdef inline bint _int32_column_append(Int32Column* col, int32_t value) nogil:
    if not _ensure_capacity(<void**>&col.data, &col.capacity, col.size + 1, sizeof(int32_t)):
        return False
    col.data[col.size] = value
    col.size += 1
    return True


cdef inline bint _uint16_column_append(UInt16Column* col, uint16_t value) nogil:
    if not _ensure_capacity(<void**>&col.data, &col.capacity, col.size + 1, sizeof(uint16_t)):
        return False
    col.data[col.size] = value
    col.size += 1
    return True


cdef inline bint _uint8_column_append(UInt8Column* col, uint8_t value) nogil:
    if not _ensure_capacity(<void**>&col.data, &col.capacity, col.size + 1, sizeof(uint8_t)):
        return False
    col.data[col.size] = value
    col.size += 1
    return True


cdef inline bint _float_column_append(FloatColumn* col, float value) nogil:
    if not _ensure_capacity(<void**>&col.data, &col.capacity, col.size + 1, sizeof(float)):
        return False
    col.data[col.size] = value
    col.size += 1
    return True


cdef inline bint _varlen_column_ensure_offsets(VarLenColumn* col, size_t needed) nogil:
    return _ensure_capacity(<void**>&col.offsets, &col.offsets_capacity, needed, sizeof(int32_t))


cdef inline bint _varlen_column_ensure_data(VarLenColumn* col, size_t needed) nogil:
    return _ensure_capacity(<void**>&col.data, &col.data_capacity, needed, sizeof(char))


cdef inline bint _varlen_column_append(VarLenColumn* col, const char* data, int32_t length) nogil:
    if col.offsets_size == 0:
        if not _varlen_column_ensure_offsets(col, 1):
            return False
        col.offsets[0] = 0
        col.offsets_size = 1

    if not _varlen_column_ensure_offsets(col, col.offsets_size + 1):
        return False
    if length < 0:
        length = 0
    if not _varlen_column_ensure_data(col, col.data_size + <size_t>length):
        return False

    if length > 0:
        memcpy(col.data + col.data_size, data, length)
        col.data_size += <size_t>length

    col.offsets[col.offsets_size] = <int32_t>col.data_size
    col.offsets_size += 1
    return True


cdef inline void _uint32_column_reset(UInt32Column* col) nogil:
    col.size = 0


cdef inline void _int32_column_reset(Int32Column* col) nogil:
    col.size = 0


cdef inline void _uint16_column_reset(UInt16Column* col) nogil:
    col.size = 0


cdef inline void _uint8_column_reset(UInt8Column* col) nogil:
    col.size = 0


cdef inline void _float_column_reset(FloatColumn* col) nogil:
    col.size = 0


cdef inline void _varlen_column_reset(VarLenColumn* col) nogil:
    col.data_size = 0
    if col.offsets_capacity == 0:
        if _varlen_column_ensure_offsets(col, 1):
            col.offsets[0] = 0
            col.offsets_size = 1
        else:
            col.offsets_size = 0
    else:
        col.offsets[0] = 0
        col.offsets_size = 1

# CIGAR operation codes from HTSlib
cdef extern from "htslib/sam.h" nogil:
    int BAM_CMATCH
    int BAM_CINS
    int BAM_CDEL
    int BAM_CREF_SKIP
    int BAM_CSOFT_CLIP
    int BAM_CHARD_CLIP
    int BAM_CPAD
    int BAM_CEQUAL
    int BAM_CDIFF
    int BAM_FUNMAP
    
    uint32_t bam_cigar_op(uint32_t c) nogil
    uint32_t bam_cigar_oplen(uint32_t c) nogil
    char bam_cigar_opchr(uint32_t c) nogil
    uint8_t* bam_get_aux(bam1_t* b)


# Compact alignment record for Parquet output
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


cdef inline void _log_progress(uint64_t processed, uint64_t written) nogil:
    if bf_should_log(1):
        bf_nogil_logf_notime(
            C_LOG_TAG,
            "processed %llu alignments (written %llu)",
            <unsigned long long>processed,
            <unsigned long long>written,
        )


cdef char* extract_cigar_string(bam1_t* b) nogil:
    """Extract CIGAR string from BAM record.
    
    Returns
    -------
    char*
        Null-terminated CIGAR string (caller must free), or NULL on error
    """
    cdef uint32_t n_cigar = b.core.n_cigar
    if n_cigar == 0:
        return <char*>NULL
    
    cdef uint32_t* cigar = bam_get_cigar(b)
    cdef char* result = <char*>malloc(<size_t>(n_cigar * 12 + 1))  # Max 11 chars per op + null
    if not result:
        return <char*>NULL
    
    cdef int pos = 0
    cdef uint32_t i, oplen
    cdef char op_chr
    cdef char* write_pos
    
    for i in range(n_cigar):
        oplen = bam_cigar_oplen(cigar[i])
        op_chr = bam_cigar_opchr(cigar[i])
        
        # Format: length + operation (e.g., "50M", "2I")
        write_pos = result + pos
        pos += snprintf(write_pos, 12, "%u%c", oplen, op_chr)
    
    result[pos] = b'\0'
    return result


cdef char* extract_sequence_string(bam1_t* b) nogil:
    """Extract sequence string from BAM record.
    
    Returns
    -------
    char*
        Null-terminated sequence string (caller must free), or NULL if empty
    """
    cdef int32_t l_seq = b.core.l_qseq
    if l_seq == 0:
        return <char*>NULL
    
    cdef uint8_t* seq = bam_get_seq(b)
    cdef char* result = <char*>malloc(<size_t>(l_seq + 1))
    if not result:
        return <char*>NULL
    
    # Nucleotide lookup table (HTSlib's standard)
    cdef char nt16_table[16]
    nt16_table[:] = [b'=', b'A', b'C', b'M', b'G', b'R', b'S', b'V', 
                      b'T', b'W', b'Y', b'H', b'K', b'D', b'B', b'N']
    cdef int i
    
    for i in range(l_seq):
        result[i] = nt16_table[bam_seqi_wrapper(seq, i)]
    
    result[l_seq] = b'\0'
    return result


cdef int populate_parquet_alignment(bam1_t* b, 
                                    sam_hdr_t* header,
                                    AlignmentScoringConfig* config,
                                    uint32_t read_id,
                                    ParquetAlignment* aln,
                                    bint include_read_name,
                                    bint include_sequence) except -1 nogil:
    """Populate ParquetAlignment structure from BAM record.
    
    Parameters
    ----------
    b : bam1_t*
        BAM record
    header : sam_hdr_t*
        BAM header
    config : AlignmentScoringConfig*
        Scoring configuration
    read_id : uint32_t
        Sequential read ID
    aln : ParquetAlignment*
        Output alignment structure
    include_read_name : bint
        Whether to extract read name (uses more space)
    include_sequence : bint
        Whether to extract sequence (uses more space)
    
    Returns
    -------
    int
        0 on success, -1 on error
    """
    cdef float pmd_score = 0.0
    cdef float* pmd_ptr = NULL
    cdef uint8_t* nm_tag
    cdef char* qname
    cdef int qname_len
    cdef uint8_t* qual
    cdef uint8_t* aux_ptr
    cdef size_t aux_len
    
    if config.calculate_pmd:
        pmd_ptr = &pmd_score
    
    # Core fields
    aln.read_id = read_id
    aln.ref_id = <uint32_t>b.core.tid
    aln.position = b.core.pos
    aln.end_position = bam_endpos(b)
    aln.mapq = b.core.qual
    aln.flag = b.core.flag
    
    # Computed metrics
    cdef double md_score = calculate_md_quality_score(b, header, config, pmd_ptr)
    aln.alignment_score = <float>md_score
    aln.pmd_score = pmd_score
    
    cdef double ani_value = 0.0
    cdef int32_t nm_value = -1
    alignment_passes_quality_filters_with_ani(b, config, &ani_value, &nm_value)
    aln.ani = <float>ani_value
    if nm_value >= 0:
        aln.num_mismatches = <uint32_t>nm_value
    else:
        aln.num_mismatches = <uint32_t>0
    
    aln.alignment_length = <uint32_t>b.core.l_qseq
    
    # Pairing information
    aln.template_length = b.core.isize
    aln.mate_ref_id = b.core.mtid
    aln.mate_position = b.core.mpos
    
    # Variable-length fields
    if include_read_name:
        qname = bam_get_qname(b)
        qname_len = b.core.l_qname
        aln.read_name = <char*>malloc(<size_t>qname_len)
        if not aln.read_name:
            return -1
        memcpy(<void*>aln.read_name, <void*>qname, <size_t>(qname_len - 1))  # Exclude null terminator
        (<char*>aln.read_name)[qname_len - 1] = 0
    else:
        aln.read_name = <char*>NULL
    
    # CIGAR string
    aln.cigar = extract_cigar_string(b)
    if b.core.n_cigar > 0 and aln.cigar == <char*>NULL:
        if aln.read_name:
            free(<void*>aln.read_name)
        return -1
    
    # Sequence (optional, uses significant space)
    if include_sequence:
        aln.sequence = extract_sequence_string(b)
        if b.core.l_qseq > 0 and aln.sequence == <char*>NULL:
            if aln.read_name:
                free(<void*>aln.read_name)
            if aln.cigar:
                free(<void*>aln.cigar)
            return -1
    else:
        aln.sequence = <char*>NULL
    
    # Quality scores (always include as BLOB)
    if b.core.l_qseq > 0:
        qual = bam_get_qual(b)
        aln.quality_length = <uint32_t>b.core.l_qseq
        aln.quality = <uint8_t*>malloc(<size_t>aln.quality_length)
        if not aln.quality:
            if aln.read_name:
                free(<void*>aln.read_name)
            if aln.cigar:
                free(<void*>aln.cigar)
            if aln.sequence:
                free(<void*>aln.sequence)
            return -1
        memcpy(<void*>aln.quality, <void*>qual, <size_t>aln.quality_length)
    else:
        aln.quality = <uint8_t*>NULL
        aln.quality_length = <uint32_t>0

    aln.tags = <uint8_t*>NULL
    aln.tags_length = 0
    aux_ptr = bam_get_aux(b)
    if aux_ptr != NULL:
        aux_len = <size_t>b.l_data - <size_t>(aux_ptr - b.data)
        if aux_len > 0:
            aln.tags = <uint8_t*>malloc(<size_t>aux_len)
            if not aln.tags:
                if aln.read_name:
                    free(<void*>aln.read_name)
                if aln.cigar:
                    free(<void*>aln.cigar)
                if aln.sequence:
                    free(<void*>aln.sequence)
                if aln.quality:
                    free(<void*>aln.quality)
                return -1
            memcpy(<void*>aln.tags, <void*>aux_ptr, <size_t>aux_len)
            aln.tags_length = <uint32_t>aux_len
    
    return 0


cdef void free_parquet_alignment(ParquetAlignment* aln) noexcept nogil:
    """Free memory allocated for ParquetAlignment."""
    if not aln:
        return
    if aln.read_name:
        free(aln.read_name)
        aln.read_name = <char*>NULL
    if aln.cigar:
        free(aln.cigar)
        aln.cigar = <char*>NULL
    if aln.sequence:
        free(aln.sequence)
        aln.sequence = <char*>NULL
    if aln.quality:
        free(aln.quality)
        aln.quality = <uint8_t*>NULL
    aln.quality_length = 0
    if aln.tags:
        free(aln.tags)
        aln.tags = <uint8_t*>NULL
    aln.tags_length = 0


cdef int process_parquet_batch(
    samFile* bam_file,
    sam_hdr_t* header,
    hts_idx_t* index,
    int64_t* reference_ids,
    int64_t ref_start,
    int64_t ref_end,
    ParquetThreadState* state,
    PartitionFlushContext* ctx,
    AlignmentScoringConfig* scoring_config,
    uint64_t* global_read_index,
    SharedCounters* counters,
    int num_partitions,
    int batch_size
) nogil:
    cdef hts_itr_t* iterator = NULL
    cdef bam1_t* bam_record = NULL
    cdef int64_t ref_idx
    cdef int64_t reference_id
    cdef int32_t ret_code
    cdef ParquetAlignment parquet_aln
    cdef int partition_id
    cdef PartitionBuffer* buffer
    cdef uint64_t processed_now = 0
    cdef uint64_t written_now = 0
    cdef uint64_t local_processed = 0
    cdef uint64_t read_block_size = 65536
    cdef uint64_t local_read_remaining = 0
    cdef uint64_t local_read_base = 0
    cdef bint include_read_names
    cdef bint include_sequences
    cdef bint include_sequence_text
    cdef bint include_pmd
    cdef int res = 0

    include_read_names = state.include_read_names
    include_sequences = state.include_sequences
    include_sequence_text = state.include_sequence_text
    include_pmd = state.include_pmd

    bam_record = bam_init1()
    if not bam_record:
        return -1

    for ref_idx in range(ref_start, ref_end):
        reference_id = reference_ids[ref_idx]
        iterator = sam_itr_queryi(index, reference_id, 0, 0x7fffffff)
        if not iterator:
            continue

        while True:
            ret_code = sam_itr_next(bam_file, iterator, bam_record)
            if ret_code < 0:
                break

            if bam_record.core.tid < 0:
                continue

            memset(&parquet_aln, 0, sizeof(ParquetAlignment))

            if global_read_index != NULL:
                if local_read_remaining == 0:
                    with gil:
                        local_read_base = global_read_index[0]
                        global_read_index[0] += read_block_size
                    local_read_remaining = read_block_size
                parquet_aln.read_id = <uint32_t>local_read_base
                local_read_base += 1
                local_read_remaining -= 1
            else:
                parquet_aln.read_id = 0

            if populate_parquet_alignment(
                bam_record,
                header,
                scoring_config,
                parquet_aln.read_id,
                &parquet_aln,
                include_read_names,
                include_sequences,
            ) != 0:
                free_parquet_alignment(&parquet_aln)
                res = -1
                break

            partition_id = (<int>parquet_aln.ref_id) % num_partitions if num_partitions > 0 else 0
            if partition_id < 0:
                partition_id = 0

            if state.active[partition_id] == 0:
                partition_buffer_initialize(
                    state.buffers + partition_id,
                    include_read_names,
                    include_sequences,
                    state.include_sequence_text,
                    include_pmd,
                )
                state.active[partition_id] = 1
            buffer = state.buffers + partition_id

            if not partition_buffer_append(buffer, &parquet_aln, include_pmd):
                free_parquet_alignment(&parquet_aln)
                res = -1
                break

            free_parquet_alignment(&parquet_aln)

            local_processed += 1
            if local_processed >= 32768:
                if counters != NULL and counters.processed != NULL:
                    with gil:
                        counters.processed[0] += local_processed
                        processed_now = counters.processed[0]
                        if counters.next_log_threshold != NULL and processed_now >= counters.next_log_threshold[0]:
                            written_now = counters.written[0] if (counters.written != NULL) else 0
                            if bf_should_log(1):
                                bf_nogil_logf_notime(
                                    C_LOG_TAG,
                                    "processed %llu alignments (written %llu)",
                                    <unsigned long long>processed_now,
                                    <unsigned long long>written_now,
                                )
                            counters.next_log_threshold[0] += 1000000
                local_processed = 0

            if batch_size > 0 and partition_buffer_row_count(buffer) >= batch_size:
                with gil:
                    flush_partition_buffer_with_gil(partition_id, buffer, ctx)

        if iterator != NULL:
            with gil:
                hts_itr_destroy(iterator)
            iterator = NULL

        if res != 0:
            break

    if iterator != NULL:
        with gil:
            hts_itr_destroy(iterator)

    if bam_record != NULL:
        with gil:
            bam_destroy1(bam_record)

    if local_processed > 0 and counters != NULL and counters.processed != NULL:
        with gil:
            counters.processed[0] += local_processed
            processed_now = counters.processed[0]
            if counters.next_log_threshold != NULL and processed_now >= counters.next_log_threshold[0]:
                written_now = counters.written[0] if (counters.written != NULL) else 0
                if bf_should_log(1):
                    bf_nogil_logf_notime(
                        C_LOG_TAG,
                        "processed %llu alignments (written %llu)",
                        <unsigned long long>processed_now,
                        <unsigned long long>written_now,
                    )
                counters.next_log_threshold[0] += 1000000

    if res == 0:
        flush_thread_buffers(state, num_partitions, ctx)
    else:
        flush_thread_buffers(state, num_partitions, ctx)

    return res

cdef object _build_reference_table(sam_hdr_t* header,
                                   int num_partitions,
                                   object pa,
                                   object ref_schema):
    """Construct a PyArrow table describing BAM references."""
    cdef int32_t total_refs = sam_hdr_nref(header)
    cdef int32_t tid
    cdef const char* name_ptr
    cdef int64_t ref_len

    ref_ids = []
    ref_names = []
    ref_lengths = []
    ref_partitions = []

    for tid in range(total_refs):
        name_ptr = sam_hdr_tid2name(header, tid)
        if name_ptr != NULL:
            ref_names.append(PyUnicode_DecodeUTF8(name_ptr, strlen(name_ptr), "strict"))
        else:
            ref_names.append("")

        ref_ids.append(tid)
        ref_len = sam_hdr_tid2len(header, tid)
        if ref_len < 0:
            ref_len = 0
        ref_lengths.append(int(ref_len))
        ref_partitions.append(tid % num_partitions)

    if bf_should_log(1):
        bf_nogil_logf_notime(
            C_LOG_TAG,
            "prepared reference table (%d references)",
            total_refs,
        )

    return pa.table(
        {
            "ref_id": ref_ids,
            "ref_name": ref_names,
            "ref_length": ref_lengths,
            "ref_partition": ref_partitions,
        },
        schema=ref_schema,
    )


cpdef object create_references_table(str bam_path, int num_partitions):
    """Return a PyArrow table describing references in the BAM header."""
    if num_partitions <= 0:
        raise ValueError("num_partitions must be a positive integer")

    import pyarrow as pa

    cdef bytes bam_path_bytes = bam_path.encode("utf-8")
    cdef samFile* bam_fp = NULL
    cdef sam_hdr_t* header = NULL

    ref_schema = pa.schema(
        [
            pa.field("ref_id", pa.uint32()),
            pa.field("ref_name", pa.string()),
            pa.field("ref_length", pa.uint32()),
            pa.field("ref_partition", pa.uint16()),
        ]
    )

    try:
        bam_fp = hts_open(bam_path_bytes, b"r")
        if not bam_fp:
            raise RuntimeError(f"Failed to open BAM file: {bam_path}")

        header = sam_hdr_read(bam_fp)
        if not header:
            raise RuntimeError("Failed to read BAM header")

        return _build_reference_table(header, num_partitions, pa, ref_schema)
    finally:
        if header != NULL:
            sam_hdr_destroy(header)
        if bam_fp != NULL:
            hts_close(bam_fp)


# Python-accessible conversion function implemented entirely in Cython/HTSlib
def convert_bam_to_parquet(str bam_path,
                          str output_base_path,
                          int num_partitions=256,
                          int batch_size=100000,
                          int num_threads=1,
                          str compression='zstd',
                          int compression_level=3,
                          bint include_read_names=True,
                          bint include_sequences=True,
                          bint include_sequence_text=False,
                          bint calculate_pmd=True,
                          int min_read_length=0,
                          int max_read_length=0,
                          float min_read_ani=0.0):
    """Convert a coordinate-sorted BAM file into partitioned Parquet format."""
    from pathlib import Path
    import json
    import os
    import pyarrow as pa
    import pyarrow.parquet as pq

    output_path = Path(output_base_path)
    output_path.mkdir(parents=True, exist_ok=True)

    alignments_dir = output_path / "alignments"
    alignments_dir.mkdir(exist_ok=True)

    references_dir = output_path / "references"
    references_dir.mkdir(exist_ok=True)

    # Ensure sequences are always emitted as requested
    if not include_sequences:
        include_sequences = True
    if include_sequence_text and not include_sequences:
        include_sequence_text = False

    cdef AlignmentScoringConfig scoring_config
    scoring_config.minimum_read_identity = 0.0
    scoring_config.minimum_read_length = 0
    scoring_config.maximum_read_length = 0
    scoring_config.global_min_score = 1e30
    scoring_config.global_max_score = -1e30
    scoring_config.calculate_pmd = calculate_pmd
    scoring_config.is_single_stranded = False

    cdef bytes bam_path_bytes = bam_path.encode("utf-8")
    cdef samFile* bam_fp = NULL
    cdef sam_hdr_t* header = NULL
    cdef hts_idx_t* bam_index = NULL

    cdef int32_t total_refs = 0
    cdef int32_t idx_tid = 0
    cdef uint64_t mapped_count = 0
    cdef uint64_t unmapped_count = 0
    cdef uint64_t total_index_alignments = 0
    cdef bint auto_partitions = (num_partitions <= 0)
    cdef bint auto_batch = (batch_size <= 0)
    cdef uint64_t target_rows = 0
    cdef uint64_t avg_partition_alignments = 0
    cdef uint64_t candidate_batch = 0
    cdef double alignments_millions = 0.0
    cdef const char* header_text = NULL
    cdef object header_py = None
    cdef object header_lower = None

    cdef int64_t* reference_ids = NULL
    cdef int64_t* reference_alignment_counts = NULL
    cdef int64_t* batch_starts = NULL
    cdef int64_t* batch_ends = NULL
    cdef int batch_count = 0
    cdef int max_batches = 0
    cdef int num_threads_c = max(1, num_threads)

    cdef samFile** thread_handles = NULL
    cdef ParquetThreadState* thread_states = NULL
    cdef uint64_t global_read_index = 0
    cdef uint64_t processed_counter = 0
    cdef uint64_t written_counter = 0
    cdef uint64_t log_threshold = 1000000
    cdef uint64_t partitions_created_counter = 0
    cdef SharedCounters counters
    counters.processed = &processed_counter
    counters.written = &written_counter
    counters.next_log_threshold = &log_threshold

    cdef PartitionFlushContext flush_ctx
    memset(&flush_ctx, 0, sizeof(PartitionFlushContext))
    flush_ctx.counters = &counters
    flush_ctx.partitions_created = &partitions_created_counter

    cdef int error_flag = 0
    cdef int* error_flag_ptr = &error_flag
    cdef int thread_id_local = 0
    cdef int batch_idx = 0
    cdef object compression_arg = None
    cdef list pending_list
    cdef object pending_entry
    cdef int pending_partition
    cdef object pending_table
    cdef uint64_t pending_rows
    cdef dict partition_writers_dict
    cdef object writer
    cdef object compression_obj
    cdef object pq_module
    cdef object align_dir
    cdef object level_arg
    cdef object pending_obj

    # Build PyArrow schemas up front
    fields = [
        pa.field("read_id", pa.uint32()),
        pa.field("ref_id", pa.uint32()),
        pa.field("position", pa.int32()),
        pa.field("end_position", pa.int32()),
        pa.field("mapq", pa.uint8()),
        pa.field("flag", pa.uint16()),
        pa.field("ani", pa.float32()),
        pa.field("alignment_score", pa.float32()),
        pa.field("pmd_score", pa.float32()),
        pa.field("num_mismatches", pa.uint32()),
        pa.field("alignment_length", pa.uint32()),
        pa.field("template_length", pa.int32()),
        pa.field("mate_ref_id", pa.int32()),
        pa.field("mate_position", pa.int32()),
    ]
    if include_read_names:
        fields.append(pa.field("read_name", pa.string()))
    fields.append(pa.field("cigar", pa.string()))
    if include_sequences:
        fields.append(pa.field("sequence_packed", pa.binary()))
        fields.append(pa.field("sequence_length", pa.uint16()))
        if include_sequence_text:
            fields.append(pa.field("sequence_text", pa.string()))
    fields.append(pa.field("quality", pa.binary()))
    fields.append(pa.field("tags", pa.binary()))
    schema = pa.schema(fields)

    ref_schema = pa.schema(
        [
            pa.field("ref_id", pa.uint32()),
            pa.field("ref_name", pa.string()),
            pa.field("ref_length", pa.uint32()),
            pa.field("ref_partition", pa.uint16()),
        ]
    )

    stats = {
        "total_alignments": 0,
        "filtered_alignments": 0,
        "written_alignments": 0,
        "partitions_created": 0,
        "include_sequences": bool(include_sequences),
        "include_sequence_text": bool(include_sequence_text),
    }
    stats["auto_num_partitions"] = bool(auto_partitions)
    stats["auto_batch_size"] = bool(auto_batch)

    partition_writers = {}

    compression = compression.lower()
    if compression not in ("zstd", "snappy", "gzip", "none"):
        raise ValueError(f"Unsupported compression codec: {compression}")

    try:
        bam_fp = hts_open(bam_path_bytes, b"r")
        if not bam_fp:
            raise RuntimeError(f"Failed to open BAM file: {bam_path}")

        header = sam_hdr_read(bam_fp)
        if not header:
            raise RuntimeError("Failed to read BAM header")

        header_text = sam_hdr_str(header)
        if header_text != NULL:
            header_py = PyUnicode_DecodeUTF8(header_text, strlen(header_text), "strict")
            header_lower = header_py.lower()
            if "so:coordinate" not in header_lower:
                raise RuntimeError("Input BAM must be coordinate-sorted (header missing SO:coordinate)")
        else:
            raise RuntimeError("BAM header missing @HD line; coordinate-sorted BAM required")

        bam_index = sam_index_load(bam_fp, bam_path_bytes)
        if bam_index == NULL:
            raise RuntimeError(
                "BAM index (.bai or .csi) is required for Parquet conversion. Run 'samtools index' first."
            )

        total_refs = sam_hdr_nref(header)
        if total_refs <= 0:
            raise RuntimeError("BAM header contains no references")

        reference_ids = <int64_t*>malloc(total_refs * sizeof(int64_t))
        reference_alignment_counts = <int64_t*>malloc(total_refs * sizeof(int64_t))
        if not reference_ids or not reference_alignment_counts:
            raise MemoryError("Failed to allocate reference arrays")

        for idx_tid in range(total_refs):
            reference_ids[idx_tid] = idx_tid
            mapped_count = 0
            unmapped_count = 0
            hts_idx_get_stat(bam_index, idx_tid, &mapped_count, &unmapped_count)
            reference_alignment_counts[idx_tid] = mapped_count
            total_index_alignments += mapped_count

        if auto_partitions:
            if total_index_alignments > 0:
                target_rows = 5000000
                num_partitions = <int>((total_index_alignments + target_rows - 1) // target_rows)
                if num_partitions < 1:
                    num_partitions = 1
                elif num_partitions > 4096:
                    num_partitions = 4096
            else:
                num_partitions = 256

        if num_partitions <= 0:
            num_partitions = 1

        if auto_batch:
            if total_index_alignments > 0:
                avg_partition_alignments = (total_index_alignments + num_partitions - 1) // num_partitions
                if avg_partition_alignments == 0:
                    avg_partition_alignments = 1
                candidate_batch = avg_partition_alignments // 16
                if candidate_batch < 50000:
                    candidate_batch = 50000
                if candidate_batch > 2000000:
                    candidate_batch = 2000000
                batch_size = <int>candidate_batch
            else:
                batch_size = 100000

        if batch_size <= 0:
            batch_size = 100000

        stats["num_partitions_used"] = num_partitions
        stats["batch_size_used"] = batch_size
        if total_index_alignments > 0:
            stats["estimated_total_alignments"] = total_index_alignments

        if total_index_alignments > 0:
            alignments_millions = (<double>total_index_alignments) / 1e6
        else:
            alignments_millions = 0.0

        if auto_partitions:
            bf_nogil_logf_notime(
                C_LOG_TAG,
                b"auto-selected %d partitions (%.2fM indexed alignments)",
                num_partitions,
                alignments_millions,
            )
        if auto_batch:
            bf_nogil_logf_notime(
                C_LOG_TAG,
                b"auto-selected batch size %d",
                batch_size,
            )

        # Write reference metadata table
        ref_table = _build_reference_table(header, num_partitions, pa, ref_schema)
        ref_output = references_dir / "references.parquet"
        if compression != "none":
            compression_arg = compression
        pq.write_table(
            ref_table,
            ref_output,
            compression=compression_arg,
            compression_level=compression_level if compression == "zstd" else None,
        )
        if bf_should_log(1):
            ref_bytes = ref_output.as_posix().encode("utf-8")
            bf_nogil_logf_notime(
                C_LOG_TAG,
                "wrote reference metadata to %s",
                <const char*>ref_bytes,
            )

        # Prepare flush context Python references (INCREF for nogil usage)
        Py_INCREF(alignments_dir)
        flush_ctx.alignments_dir = <PyObject*>alignments_dir
        Py_INCREF(schema)
        flush_ctx.schema = <PyObject*>schema
        Py_INCREF(compression)
        flush_ctx.compression = <PyObject*>compression
        flush_ctx.compression_level = compression_level
        flush_ctx.batch_size_hint = batch_size
        Py_INCREF(pa)
        flush_ctx.pa_module = <PyObject*>pa
        Py_INCREF(pq)
        flush_ctx.pq_module = <PyObject*>pq
        Py_INCREF(partition_writers)
        flush_ctx.partition_writers = <PyObject*>partition_writers
        pending_obj = PyList_New(0)
        if pending_obj is None:
            raise MemoryError("Failed to allocate pending table list")
        flush_ctx.pending_tables = <PyObject*>pending_obj

        # Allocate batches
        max_batches = total_refs
        if max_batches > 1024:
            max_batches = 1024
        if max_batches < num_threads_c:
            max_batches = num_threads_c
        batch_starts = <int64_t*>malloc(max_batches * sizeof(int64_t))
        batch_ends = <int64_t*>malloc(max_batches * sizeof(int64_t))
        if not batch_starts or not batch_ends:
            raise MemoryError("Failed to allocate batch arrays")

        batch_count = create_balanced_batches_greedy(
            reference_ids,
            reference_alignment_counts,
            total_refs,
            max_batches,
            batch_starts,
            batch_ends,
            num_threads_c,
            False
        )
        if batch_count <= 0:
            raise RuntimeError("No batches created for Parquet streaming")

        # Open per-thread BAM handles
        thread_handles = <samFile**>malloc(num_threads_c * sizeof(samFile*))
        if not thread_handles:
            raise MemoryError("Failed to allocate thread handle array")
        for idx_tid in range(num_threads_c):
            thread_handles[idx_tid] = NULL
        for idx_tid in range(num_threads_c):
            thread_handles[idx_tid] = hts_open(bam_path_bytes, b"r")
            if not thread_handles[idx_tid]:
                raise RuntimeError(f"Failed to open BAM file handle for thread {idx_tid}")

        # Initialize per-thread Parquet buffers
        thread_states = <ParquetThreadState*>malloc(num_threads_c * sizeof(ParquetThreadState))
        if not thread_states:
            raise MemoryError("Failed to allocate thread state array")
        memset(thread_states, 0, num_threads_c * sizeof(ParquetThreadState))
        for idx_tid in range(num_threads_c):
            if parquet_thread_state_init(thread_states + idx_tid, num_partitions, include_read_names, include_sequences, include_sequence_text, calculate_pmd) != 0:
                raise MemoryError(f"Failed to initialize thread state {idx_tid}")

        # Parallel streaming
        with nogil:
            for batch_idx in prange(batch_count, num_threads=num_threads_c, schedule='static'):
                if error_flag_ptr[0] != 0:
                    continue
                thread_id_local = threadid()
                if process_parquet_batch(
                    thread_handles[thread_id_local],
                    header,
                    bam_index,
                    reference_ids,
                    batch_starts[batch_idx],
                    batch_ends[batch_idx],
                    thread_states + thread_id_local,
                    &flush_ctx,
                    &scoring_config,
                    &global_read_index,
                    &counters,
                    num_partitions,
                    batch_size
                ) != 0:
                    error_flag_ptr[0] = 1

        if error_flag != 0:
            raise MemoryError("Parquet streaming aborted due to allocation failure")

        pending_list = <list>flush_ctx.pending_tables
        partition_writers_dict = <dict>flush_ctx.partition_writers
        pq_module = <object>flush_ctx.pq_module
        align_dir = <object>flush_ctx.alignments_dir
        for pending_entry in pending_list:
            pending_partition = <int>pending_entry[0]
            pending_table = pending_entry[1]
            pending_rows = <uint64_t>pending_entry[2]

            writer = partition_writers_dict.get(pending_partition)
            if writer is None:
                partition_dir = align_dir / f"ref_partition={pending_partition:04d}"
                partition_dir.mkdir(exist_ok=True)
                output_file = partition_dir / "data.parquet"
                compression_obj = <object>flush_ctx.compression
                if compression_obj == "none":
                    compression_obj = None
                level_arg = None
                if compression_obj == "zstd" and flush_ctx.compression_level > 0:
                    level_arg = flush_ctx.compression_level
                writer = (<object>pq_module).ParquetWriter(
                    output_file,
                    <object>flush_ctx.schema,
                    compression=compression_obj,
                    compression_level=level_arg,
                    use_dictionary=True,
                    write_statistics=True,
                )
                partition_writers_dict[pending_partition] = writer
                if flush_ctx.partitions_created != NULL:
                    flush_ctx.partitions_created[0] += 1
                bf_nogil_logf_notime(
                    C_LOG_TAG,
                    "created partition %04d",
                    pending_partition,
                )

            if flush_ctx.batch_size_hint > 0:
                writer.write_table(pending_table, row_group_size=flush_ctx.batch_size_hint)
            else:
                writer.write_table(pending_table)

            if flush_ctx.counters != NULL and flush_ctx.counters.written != NULL:
                flush_ctx.counters.written[0] += pending_rows

        pending_list.clear()

    finally:
        # Close Parquet writers
        for writer in partition_writers.values():
            writer.close()

        # Release flush context Python references
        if flush_ctx.partition_writers != NULL:
            Py_DECREF(<object>flush_ctx.partition_writers)
            flush_ctx.partition_writers = NULL
        if flush_ctx.pq_module != NULL:
            Py_DECREF(<object>flush_ctx.pq_module)
            flush_ctx.pq_module = NULL
        if flush_ctx.pa_module != NULL:
            Py_DECREF(<object>flush_ctx.pa_module)
            flush_ctx.pa_module = NULL
        if flush_ctx.compression != NULL:
            Py_DECREF(<object>flush_ctx.compression)
            flush_ctx.compression = NULL
        if flush_ctx.schema != NULL:
            Py_DECREF(<object>flush_ctx.schema)
            flush_ctx.schema = NULL
        if flush_ctx.alignments_dir != NULL:
            Py_DECREF(<object>flush_ctx.alignments_dir)
            flush_ctx.alignments_dir = NULL
        if flush_ctx.pending_tables != NULL:
            Py_DECREF(<object>flush_ctx.pending_tables)
            flush_ctx.pending_tables = NULL

        if thread_states != NULL:
            for idx_tid in range(num_threads_c):
                parquet_thread_state_free(thread_states + idx_tid, num_partitions)
            free(thread_states)
            thread_states = NULL

        if thread_handles != NULL:
            for idx_tid in range(num_threads_c):
                if thread_handles[idx_tid]:
                    hts_close(thread_handles[idx_tid])
            free(thread_handles)
            thread_handles = NULL

        if batch_starts != NULL:
            free(batch_starts)
            batch_starts = NULL
        if batch_ends != NULL:
            free(batch_ends)
            batch_ends = NULL

        if reference_ids != NULL:
            free(reference_ids)
            reference_ids = NULL
        if reference_alignment_counts != NULL:
            free(reference_alignment_counts)
            reference_alignment_counts = NULL

        if bam_index != NULL:
            hts_idx_destroy(bam_index)
            bam_index = NULL
        if header != NULL:
            sam_hdr_destroy(header)
            header = NULL
        if bam_fp != NULL:
            hts_close(bam_fp)
            bam_fp = NULL

    if error_flag != 0:
        raise MemoryError("Parquet conversion failed")

    stats["total_alignments"] = int(processed_counter)
    stats["written_alignments"] = int(written_counter)
    stats["partitions_created"] = int(partitions_created_counter)

    metadata = {
        "source_bam": os.path.abspath(bam_path),
        "num_partitions": num_partitions,
        "batch_size": batch_size,
        "compression": compression,
        "compression_level": compression_level,
        "include_read_names": bool(include_read_names),
        "include_sequences": bool(include_sequences),
        "include_sequence_text": bool(include_sequence_text),
        "calculate_pmd": bool(calculate_pmd),
        "auto_num_partitions": bool(auto_partitions),
        "auto_batch_size": bool(auto_batch),
        "estimated_total_alignments": int(total_index_alignments) if total_index_alignments > 0 else None,
        "statistics": stats,
    }

    metadata_path = output_path / "_metadata.json"
    with open(metadata_path, "w") as fh:
        json.dump(metadata, fh, indent=2)

    return stats
cdef struct PartitionBuffer:
    UInt32Column read_id
    UInt32Column ref_id
    Int32Column position
    Int32Column end_position
    UInt8Column mapq
    UInt16Column flag
    FloatColumn ani
    FloatColumn alignment_score
    FloatColumn pmd_score
    UInt32Column num_mismatches
    UInt32Column alignment_length
    Int32Column template_length
    Int32Column mate_ref_id
    Int32Column mate_position
    VarLenColumn read_name
    VarLenColumn cigar
    VarLenColumn sequence
    UInt16Column sequence_length
    VarLenColumn sequence_text
    VarLenColumn quality
    VarLenColumn tags
    bint include_read_names
    bint include_sequences
    bint include_sequence_text
    bint include_pmd
    bint initialized


cdef inline void partition_buffer_initialize(PartitionBuffer* buf,
                                             bint include_read_names,
                                             bint include_sequences,
                                             bint include_sequence_text,
                                             bint include_pmd) nogil:
    if buf.initialized:
        return
    _uint32_column_init(&buf.read_id)
    _uint32_column_init(&buf.ref_id)
    _int32_column_init(&buf.position)
    _int32_column_init(&buf.end_position)
    _uint8_column_init(&buf.mapq)
    _uint16_column_init(&buf.flag)
    _float_column_init(&buf.ani)
    _float_column_init(&buf.alignment_score)
    _float_column_init(&buf.pmd_score)
    _uint32_column_init(&buf.num_mismatches)
    _uint32_column_init(&buf.alignment_length)
    _int32_column_init(&buf.template_length)
    _int32_column_init(&buf.mate_ref_id)
    _int32_column_init(&buf.mate_position)
    _varlen_column_init(&buf.read_name)
    _varlen_column_init(&buf.cigar)
    _varlen_column_init(&buf.sequence)
    _uint16_column_init(&buf.sequence_length)
    _varlen_column_init(&buf.sequence_text)
    _varlen_column_init(&buf.quality)
    _varlen_column_init(&buf.tags)
    buf.include_read_names = include_read_names
    buf.include_sequences = include_sequences
    buf.include_sequence_text = include_sequence_text
    buf.include_pmd = include_pmd
    buf.initialized = True


cdef inline void partition_buffer_free(PartitionBuffer* buf) nogil:
    if not buf.initialized:
        return
    _uint32_column_free(&buf.read_id)
    _uint32_column_free(&buf.ref_id)
    _int32_column_free(&buf.position)
    _int32_column_free(&buf.end_position)
    _uint8_column_free(&buf.mapq)
    _uint16_column_free(&buf.flag)
    _float_column_free(&buf.ani)
    _float_column_free(&buf.alignment_score)
    _float_column_free(&buf.pmd_score)
    _uint32_column_free(&buf.num_mismatches)
    _uint32_column_free(&buf.alignment_length)
    _int32_column_free(&buf.template_length)
    _int32_column_free(&buf.mate_ref_id)
    _int32_column_free(&buf.mate_position)
    _varlen_column_free(&buf.read_name)
    _varlen_column_free(&buf.cigar)
    _varlen_column_free(&buf.sequence)
    _uint16_column_free(&buf.sequence_length)
    _varlen_column_free(&buf.sequence_text)
    _varlen_column_free(&buf.quality)
    _varlen_column_free(&buf.tags)
    buf.initialized = False


cdef inline bint partition_buffer_append(PartitionBuffer* buf,
                                         ParquetAlignment* aln,
                                         bint calculate_pmd) nogil:
    if not buf.initialized:
        return False
    if not _uint32_column_append(&buf.read_id, aln.read_id):
        return False
    if not _uint32_column_append(&buf.ref_id, aln.ref_id):
        return False
    if not _int32_column_append(&buf.position, aln.position):
        return False
    if not _int32_column_append(&buf.end_position, aln.end_position):
        return False
    if not _uint8_column_append(&buf.mapq, aln.mapq):
        return False
    if not _uint16_column_append(&buf.flag, aln.flag):
        return False
    if not _float_column_append(&buf.ani, aln.ani):
        return False
    if not _float_column_append(&buf.alignment_score, aln.alignment_score):
        return False
    if calculate_pmd:
        if not _float_column_append(&buf.pmd_score, aln.pmd_score):
            return False
    else:
        if not _float_column_append(&buf.pmd_score, <float>NAN):
            return False
    if not _uint32_column_append(&buf.num_mismatches, aln.num_mismatches):
        return False
    if not _uint32_column_append(&buf.alignment_length, aln.alignment_length):
        return False
    if not _int32_column_append(&buf.template_length, aln.template_length):
        return False
    if not _int32_column_append(&buf.mate_ref_id, aln.mate_ref_id):
        return False
    if not _int32_column_append(&buf.mate_position, aln.mate_position):
        return False

    cdef int32_t cigar_len = 0
    cdef int32_t seq_len = 0
    cdef int32_t name_len = 0
    cdef const char* cigar_ptr = NULL
    cdef const char* seq_ptr = NULL
    cdef char* packed_ptr = NULL
    cdef uint32_t packed_len = 0

    if buf.include_read_names:
        if aln.read_name != NULL:
            name_len = <int32_t>strlen(aln.read_name)
            if not _varlen_column_append(&buf.read_name, aln.read_name, name_len):
                return False
        else:
            if not _varlen_column_append(&buf.read_name, NULL, 0):
                return False

    if aln.cigar != NULL:
        cigar_len = <int32_t>strlen(aln.cigar)
        cigar_ptr = aln.cigar
    else:
        cigar_ptr = NULL
    if not _varlen_column_append(&buf.cigar, cigar_ptr, cigar_len):
        return False

    if buf.include_sequences:
        if aln.sequence != NULL and aln.alignment_length > 0:
            if not pack_sequence_bases(aln.sequence, <int32_t>aln.alignment_length, &packed_ptr, &packed_len):
                return False
            if not _varlen_column_append(&buf.sequence, packed_ptr, packed_len):
                if packed_ptr != NULL:
                    free(packed_ptr)
                    packed_ptr = NULL
                return False
            if packed_ptr != NULL:
                free(packed_ptr)
                packed_ptr = NULL
        else:
            if not _varlen_column_append(&buf.sequence, NULL, 0):
                return False
        if not _uint16_column_append(&buf.sequence_length, <uint16_t>aln.alignment_length):
            return False
        if buf.include_sequence_text:
            if aln.sequence != NULL and aln.alignment_length > 0:
                seq_len = <int32_t>strlen(aln.sequence)
                if not _varlen_column_append(&buf.sequence_text, aln.sequence, seq_len):
                    return False
            else:
                if not _varlen_column_append(&buf.sequence_text, NULL, 0):
                    return False

    if aln.quality != NULL and aln.quality_length > 0:
        if not _varlen_column_append(&buf.quality, <const char*>aln.quality, aln.quality_length):
            return False
    else:
        if not _varlen_column_append(&buf.quality, NULL, 0):
            return False

    if aln.tags != NULL and aln.tags_length > 0:
        if not _varlen_column_append(&buf.tags, <const char*>aln.tags, <int32_t>aln.tags_length):
            return False
    else:
        if not _varlen_column_append(&buf.tags, NULL, 0):
            return False

    return True


cdef inline size_t partition_buffer_row_count(PartitionBuffer* buf) nogil:
    if not buf.initialized:
        return 0
    return buf.read_id.size


cdef inline void partition_buffer_reset(PartitionBuffer* buf) nogil:
    if not buf.initialized:
        return
    _uint32_column_reset(&buf.read_id)
    _uint32_column_reset(&buf.ref_id)
    _int32_column_reset(&buf.position)
    _int32_column_reset(&buf.end_position)
    _uint8_column_reset(&buf.mapq)
    _uint16_column_reset(&buf.flag)
    _float_column_reset(&buf.ani)
    _float_column_reset(&buf.alignment_score)
    _float_column_reset(&buf.pmd_score)
    _uint32_column_reset(&buf.num_mismatches)
    _uint32_column_reset(&buf.alignment_length)
    _int32_column_reset(&buf.template_length)
    _int32_column_reset(&buf.mate_ref_id)
    _int32_column_reset(&buf.mate_position)
    if buf.include_read_names:
        _varlen_column_reset(&buf.read_name)
    _varlen_column_reset(&buf.cigar)
    if buf.include_sequences:
        _varlen_column_reset(&buf.sequence)
        _uint16_column_reset(&buf.sequence_length)
        if buf.include_sequence_text:
            _varlen_column_reset(&buf.sequence_text)
    _varlen_column_reset(&buf.quality)
    _varlen_column_reset(&buf.tags)
cdef inline object _buffer_from_ptr(object pa, void* ptr, Py_ssize_t size):
    if size <= 0:
        size = 0
    cdef object py_bytes = PyBytes_FromStringAndSize(<char*>ptr, size)
    if py_bytes is None:
        raise MemoryError()
    return pa.py_buffer(py_bytes)


cdef object partition_buffer_to_arrow_table(PartitionBuffer* buf,
                                            object pa,
                                            object schema):
    cdef size_t n
    cdef list arrays
    cdef object arr
    cdef object offsets_buf
    cdef object data_buf

    if buf is NULL or not buf.initialized:
        return None

    n = buf.read_id.size
    if n == 0:
        return None

    arrays = []

    data_buf = _buffer_from_ptr(pa, <void*>buf.read_id.data, <Py_ssize_t>(n * sizeof(uint32_t)))
    arrays.append(pa.Array.from_buffers(pa.uint32(), n, [None, data_buf]))

    data_buf = _buffer_from_ptr(pa, <void*>buf.ref_id.data, <Py_ssize_t>(n * sizeof(uint32_t)))
    arrays.append(pa.Array.from_buffers(pa.uint32(), n, [None, data_buf]))

    data_buf = _buffer_from_ptr(pa, <void*>buf.position.data, <Py_ssize_t>(n * sizeof(int32_t)))
    arrays.append(pa.Array.from_buffers(pa.int32(), n, [None, data_buf]))

    data_buf = _buffer_from_ptr(pa, <void*>buf.end_position.data, <Py_ssize_t>(n * sizeof(int32_t)))
    arrays.append(pa.Array.from_buffers(pa.int32(), n, [None, data_buf]))

    data_buf = _buffer_from_ptr(pa, <void*>buf.mapq.data, <Py_ssize_t>(n * sizeof(uint8_t)))
    arrays.append(pa.Array.from_buffers(pa.uint8(), n, [None, data_buf]))

    data_buf = _buffer_from_ptr(pa, <void*>buf.flag.data, <Py_ssize_t>(n * sizeof(uint16_t)))
    arrays.append(pa.Array.from_buffers(pa.uint16(), n, [None, data_buf]))

    data_buf = _buffer_from_ptr(pa, <void*>buf.ani.data, <Py_ssize_t>(n * sizeof(float)))
    arrays.append(pa.Array.from_buffers(pa.float32(), n, [None, data_buf]))

    data_buf = _buffer_from_ptr(pa, <void*>buf.alignment_score.data, <Py_ssize_t>(n * sizeof(float)))
    arrays.append(pa.Array.from_buffers(pa.float32(), n, [None, data_buf]))

    data_buf = _buffer_from_ptr(pa, <void*>buf.pmd_score.data, <Py_ssize_t>(n * sizeof(float)))
    arrays.append(pa.Array.from_buffers(pa.float32(), n, [None, data_buf]))

    data_buf = _buffer_from_ptr(pa, <void*>buf.num_mismatches.data, <Py_ssize_t>(n * sizeof(uint32_t)))
    arrays.append(pa.Array.from_buffers(pa.uint32(), n, [None, data_buf]))

    data_buf = _buffer_from_ptr(pa, <void*>buf.alignment_length.data, <Py_ssize_t>(n * sizeof(uint32_t)))
    arrays.append(pa.Array.from_buffers(pa.uint32(), n, [None, data_buf]))

    data_buf = _buffer_from_ptr(pa, <void*>buf.template_length.data, <Py_ssize_t>(n * sizeof(int32_t)))
    arrays.append(pa.Array.from_buffers(pa.int32(), n, [None, data_buf]))

    data_buf = _buffer_from_ptr(pa, <void*>buf.mate_ref_id.data, <Py_ssize_t>(n * sizeof(int32_t)))
    arrays.append(pa.Array.from_buffers(pa.int32(), n, [None, data_buf]))

    data_buf = _buffer_from_ptr(pa, <void*>buf.mate_position.data, <Py_ssize_t>(n * sizeof(int32_t)))
    arrays.append(pa.Array.from_buffers(pa.int32(), n, [None, data_buf]))

    if buf.include_read_names:
        if buf.read_name.offsets_size != n + 1:
            raise ValueError("read_name offsets length mismatch")
        offsets_buf = _buffer_from_ptr(pa, <void*>buf.read_name.offsets, <Py_ssize_t>(buf.read_name.offsets_size * sizeof(int32_t)))
        data_buf = _buffer_from_ptr(pa, <void*>buf.read_name.data, <Py_ssize_t>buf.read_name.data_size)
        arr = pa.Array.from_buffers(pa.string(), n, [None, offsets_buf, data_buf])
        arrays.append(arr)

    if buf.cigar.offsets_size != n + 1:
        raise ValueError("cigar offsets length mismatch")
    offsets_buf = _buffer_from_ptr(pa, <void*>buf.cigar.offsets, <Py_ssize_t>(buf.cigar.offsets_size * sizeof(int32_t)))
    data_buf = _buffer_from_ptr(pa, <void*>buf.cigar.data, <Py_ssize_t>buf.cigar.data_size)
    arr = pa.Array.from_buffers(pa.string(), n, [None, offsets_buf, data_buf])
    arrays.append(arr)

    if buf.include_sequences:
        if buf.sequence.offsets_size != n + 1:
            raise ValueError("sequence offsets length mismatch")
        offsets_buf = _buffer_from_ptr(pa, <void*>buf.sequence.offsets, <Py_ssize_t>(buf.sequence.offsets_size * sizeof(int32_t)))
        data_buf = _buffer_from_ptr(pa, <void*>buf.sequence.data, <Py_ssize_t>buf.sequence.data_size)
        arrays.append(pa.Array.from_buffers(pa.binary(), n, [None, offsets_buf, data_buf]))

        data_buf = _buffer_from_ptr(pa, <void*>buf.sequence_length.data, <Py_ssize_t>(n * sizeof(uint16_t)))
        arrays.append(pa.Array.from_buffers(pa.uint16(), n, [None, data_buf]))

        if buf.include_sequence_text:
            if buf.sequence_text.offsets_size != n + 1:
                raise ValueError("sequence_text offsets length mismatch")
            offsets_buf = _buffer_from_ptr(pa, <void*>buf.sequence_text.offsets, <Py_ssize_t>(buf.sequence_text.offsets_size * sizeof(int32_t)))
            data_buf = _buffer_from_ptr(pa, <void*>buf.sequence_text.data, <Py_ssize_t>buf.sequence_text.data_size)
            arr = pa.Array.from_buffers(pa.string(), n, [None, offsets_buf, data_buf])
            arrays.append(arr)

    if buf.quality.offsets_size != n + 1:
        raise ValueError("quality offsets length mismatch")
    offsets_buf = _buffer_from_ptr(pa, <void*>buf.quality.offsets, <Py_ssize_t>(buf.quality.offsets_size * sizeof(int32_t)))
    data_buf = _buffer_from_ptr(pa, <void*>buf.quality.data, <Py_ssize_t>buf.quality.data_size)
    arr = pa.Array.from_buffers(pa.binary(), n, [None, offsets_buf, data_buf])
    arrays.append(arr)

    if buf.tags.offsets_size != n + 1:
        raise ValueError("tags offsets length mismatch")
    offsets_buf = _buffer_from_ptr(pa, <void*>buf.tags.offsets, <Py_ssize_t>(buf.tags.offsets_size * sizeof(int32_t)))
    data_buf = _buffer_from_ptr(pa, <void*>buf.tags.data, <Py_ssize_t>buf.tags.data_size)
    arr = pa.Array.from_buffers(pa.binary(), n, [None, offsets_buf, data_buf])
    arrays.append(arr)

    return pa.Table.from_arrays(arrays, schema=schema)


cdef struct SharedCounters:
    uint64_t* processed
    uint64_t* written
    uint64_t* next_log_threshold


cdef struct PartitionFlushContext:
    PyObject* alignments_dir
    PyObject* schema
    PyObject* compression
    int compression_level
    int batch_size_hint
    PyObject* pa_module
    PyObject* pq_module
    PyObject* partition_writers
    SharedCounters* counters
    uint64_t* partitions_created
    PyObject* pending_tables


cdef struct ParquetThreadState:
    PartitionBuffer* buffers
    uint8_t* active
    bint include_read_names
    bint include_sequences
    bint include_sequence_text
    bint include_pmd


cdef int parquet_thread_state_init(ParquetThreadState* state,
                                   int num_partitions,
                                   bint include_read_names,
                                   bint include_sequences,
                                   bint include_sequence_text,
                                   bint include_pmd) nogil:
    state.buffers = <PartitionBuffer*>malloc(num_partitions * sizeof(PartitionBuffer))
    if not state.buffers:
        return -1
    state.active = <uint8_t*>calloc(num_partitions, sizeof(uint8_t))
    if not state.active:
        free(state.buffers)
        state.buffers = NULL
        return -1
    state.include_read_names = include_read_names
    state.include_sequences = include_sequences
    state.include_sequence_text = include_sequence_text
    state.include_pmd = include_pmd
    cdef int i
    for i in range(num_partitions):
        memset(state.buffers + i, 0, sizeof(PartitionBuffer))
    return 0


cdef void parquet_thread_state_free(ParquetThreadState* state,
                                    int num_partitions) nogil:
    cdef int i
    if state.buffers != NULL:
        for i in range(num_partitions):
            if state.active != NULL and state.active[i]:
                partition_buffer_free(state.buffers + i)
        free(state.buffers)
        state.buffers = NULL
    if state.active != NULL:
        free(state.active)
        state.active = NULL


cdef int flush_thread_buffers(ParquetThreadState* state,
                              int num_partitions,
                              PartitionFlushContext* ctx) nogil:
    cdef int i
    cdef uint64_t _rows
    if state is NULL or ctx is NULL:
        return 0
    for i in range(num_partitions):
        if state.active != NULL and state.active[i]:
            with gil:
                _rows = flush_partition_buffer_with_gil(i, state.buffers + i, ctx)
            # rows accounted inside flush
    return 0


cdef uint64_t flush_partition_buffer_with_gil(int partition_id,
                                              PartitionBuffer* buffer,
                                              PartitionFlushContext* ctx) with gil:
    if buffer is NULL or not buffer.initialized:
        return 0

    cdef size_t row_count = buffer.read_id.size
    cdef object pa = <object>ctx.pa_module
    cdef object schema = <object>ctx.schema

    bf_nogil_logf_notime(
        C_LOG_TAG,
        "flush partition=%d rows=%zu name_offsets=%zu cigar_offsets=%zu seq_offsets=%zu seq_text_offsets=%zu q_offsets=%zu tag_offsets=%zu",
        partition_id,
        row_count,
        buffer.read_name.offsets_size if buffer.include_read_names else 0,
        buffer.cigar.offsets_size,
        buffer.sequence.offsets_size if buffer.include_sequences else 0,
        buffer.sequence_text.offsets_size if buffer.include_sequence_text else 0,
        buffer.quality.offsets_size,
        buffer.tags.offsets_size,
    )

    cdef object table = partition_buffer_to_arrow_table(buffer, pa, schema)
    if table is None:
        return 0

    table = table.combine_chunks()
    table.validate(full=True)

    cdef uint64_t rows = <uint64_t>table.num_rows

    cdef list pending = <list>ctx.pending_tables
    pending.append((partition_id, table, rows))

    partition_buffer_reset(buffer)
    return rows
