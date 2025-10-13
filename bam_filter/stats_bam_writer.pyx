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

from libc.stdint cimport int32_t, int64_t, uint32_t, uint8_t, uint16_t
from libc.stdlib cimport malloc, free

from bam_filter.stats cimport RefStats, FilterConditions
from bam_filter.processor cimport (
    min_int32,
    samFile,
    sam_hdr_t,
)
from bam_filter.processor_mapping cimport ReferenceMapping, create_filtered_header_efficient

from bam_filter.stats_bam_writer cimport ReferenceFilter

cdef extern from "time.h":
    cdef struct timespec:
        long tv_sec
        long tv_nsec
    int clock_gettime(int clk_id, timespec* tp) nogil
    int CLOCK_MONOTONIC

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil

cdef extern from "htslib/sam.h":
    int bam_aux2i(const uint8_t *s) nogil
    ctypedef struct bam1_core_t:
        int32_t tid
        int64_t pos
        uint8_t qual
        uint16_t flag
        uint16_t l_qname
        uint32_t n_cigar
        int32_t l_qseq
    ctypedef struct bam1_t:
        bam1_core_t core
        uint8_t* data
        int l_data
        uint32_t m_data
    ctypedef struct hts_idx_t
    ctypedef struct hts_itr_t

    samFile* hts_open(const char* fn, const char* mode) nogil
    int hts_close(samFile* fp) nogil
    int sam_hdr_write(samFile* fp, const sam_hdr_t* h) nogil
    int sam_write1(samFile* fp, const sam_hdr_t* h, const bam1_t* b) nogil
    hts_idx_t* sam_index_load(samFile* fp, const char* fn) nogil
    hts_itr_t* sam_itr_queryi(hts_idx_t* idx, int tid, int beg, int end) nogil
    int sam_itr_next(samFile* fp, hts_itr_t* iter, bam1_t* b) nogil
    void hts_itr_destroy(hts_itr_t* iter) nogil
    void hts_idx_destroy(hts_idx_t* idx) nogil
    bam1_t* bam_init1() nogil
    void bam_destroy1(bam1_t* b) nogil
    uint8_t* bam_aux_get(bam1_t* b, const char* tag) nogil
    void sam_hdr_destroy(sam_hdr_t* h) nogil
    int hts_set_threads(samFile* fp, int n) nogil

cdef extern from "htslib/hts.h":
    int hts_set_opt(void* fp, int opt, ...) nogil
    cdef int HTS_OPT_CACHE_SIZE


cdef bint passes_filters(RefStats* stats, FilterConditions* filters) noexcept nogil:
    """Check if reference passes all ENABLED filter conditions."""
    cdef double cov_evenness_tmp

    if filters.enable_min_coverage_mean and stats.mean_coverage < filters.min_coverage_mean:
        return False

    if filters.enable_min_avg_read_ani and stats.ani_mean < filters.min_avg_read_ani:
        return False

    if filters.enable_min_breadth and stats.breadth < filters.min_breadth:
        return False

    if filters.enable_min_expected_breadth_ratio and stats.breadth_exp_ratio < filters.min_expected_breadth_ratio:
        return False

    if filters.enable_min_coverage_evenness:
        cov_evenness_tmp = stats.cov_evenness if stats.mean_coverage >= 1.0 else 1.0
        if cov_evenness_tmp < filters.min_coverage_evenness:
            return False

    if filters.enable_max_coeff_var and stats.c_v > filters.max_coeff_var:
        return False

    if filters.enable_min_norm_entropy and stats.norm_entropy < filters.min_norm_entropy:
        return False

    if filters.enable_max_norm_gini and stats.norm_gini > filters.max_norm_gini:
        return False

    return True


cdef ReferenceFilter* create_reference_filter(
    RefStats* global_ref_stats,
    FilterConditions* filters,
    int n_refs
) noexcept nogil:
    cdef ReferenceFilter* ref_filter = <ReferenceFilter*>malloc(sizeof(ReferenceFilter))
    if ref_filter == NULL:
        return NULL

    ref_filter.n_total_refs = n_refs
    ref_filter.tid_mapping = <int32_t*>malloc(n_refs * sizeof(int32_t))
    if ref_filter.tid_mapping == NULL:
        free(ref_filter)
        return NULL

    cdef int i
    for i in range(n_refs):
        ref_filter.tid_mapping[i] = -1

    cdef int32_t new_tid = 0
    cdef RefStats* stats

    for i in range(n_refs):
        stats = &global_ref_stats[i]
        if stats.n_alns > 0 and passes_filters(stats, filters):
            ref_filter.tid_mapping[i] = new_tid
            new_tid += 1

    ref_filter.n_filtered_refs = new_tid

    if new_tid > 0:
        ref_filter.reverse_mapping = <int32_t*>malloc(new_tid * sizeof(int32_t))
        if ref_filter.reverse_mapping == NULL:
            free(ref_filter.tid_mapping)
            free(ref_filter)
            return NULL

        for i in range(n_refs):
            if ref_filter.tid_mapping[i] >= 0:
                ref_filter.reverse_mapping[ref_filter.tid_mapping[i]] = i
    else:
        ref_filter.reverse_mapping = NULL

    return ref_filter


cdef void destroy_reference_filter(ReferenceFilter* ref_filter) noexcept nogil:
    if ref_filter != NULL:
        if ref_filter.tid_mapping != NULL:
            free(ref_filter.tid_mapping)
        if ref_filter.reverse_mapping != NULL:
            free(ref_filter.reverse_mapping)
        free(ref_filter)


cdef int write_filtered_bam_streaming(
    samFile* input_bam,
    hts_idx_t* existing_index,
    const char* input_bam_path,
    const char* output_bam_path,
    sam_hdr_t* original_header,
    ReferenceFilter* ref_filter,
    int num_threads,
    double min_read_ani_c,
    int min_read_length_c,
    int max_read_length_c
) except -1 nogil:
    cdef samFile* output_bam = NULL
    cdef hts_idx_t* bam_index = NULL
    cdef sam_hdr_t* filtered_header = NULL
    cdef ReferenceMapping* mapping = NULL
    cdef bam1_t* record = NULL
    cdef hts_itr_t* iterator = NULL
    cdef int ret = 0
    cdef int64_t alignments_written = 0
    cdef int64_t alignments_processed = 0
    cdef int32_t old_tid, new_tid
    cdef int32_t orig_tid
    cdef int64_t ref_alignments = 0
    cdef uint32_t mapped_new
    cdef int cache_size = 256 * 1024 * 1024
    cdef int i
    cdef int32_t aln_read_length
    cdef uint8_t* nm_tag_local
    cdef int nm_val_local
    cdef double aln_ani
    cdef timespec ts_write_start, ts_write_end
    cdef double write_sec = 0.0

    if ref_filter.n_filtered_refs == 0:
        bf_nogil_logf_notime(NULL, "No references pass filtering - no BAM output created\n")
        return 0

    mapping = <ReferenceMapping*>malloc(sizeof(ReferenceMapping))
    if mapping == NULL:
        return -1

    mapping.n_retained_refs = ref_filter.n_filtered_refs
    mapping.n_original_refs = ref_filter.n_total_refs

    mapping.new_to_old_tid = <uint32_t*>malloc(mapping.n_retained_refs * sizeof(uint32_t))
    mapping.old_to_new_tid = <uint32_t*>malloc(mapping.n_original_refs * sizeof(uint32_t))

    if mapping.new_to_old_tid == NULL or mapping.old_to_new_tid == NULL:
        if mapping.new_to_old_tid != NULL:
            free(mapping.new_to_old_tid)
        if mapping.old_to_new_tid != NULL:
            free(mapping.old_to_new_tid)
        free(mapping)
        return -1

    for i in range(mapping.n_original_refs):
        mapping.old_to_new_tid[i] = <uint32_t>(-1)

    for new_tid in range(mapping.n_retained_refs):
        old_tid = ref_filter.reverse_mapping[new_tid]
        mapping.new_to_old_tid[new_tid] = <uint32_t>old_tid
        mapping.old_to_new_tid[<uint32_t>old_tid] = <uint32_t>new_tid

    filtered_header = create_filtered_header_efficient(original_header, mapping)

    if filtered_header == NULL:
        free(mapping.new_to_old_tid)
        free(mapping.old_to_new_tid)
        free(mapping)
        return -1

    cdef bint opened_input_bam_here = False

    if input_bam == NULL:
        input_bam = hts_open(input_bam_path, b"r")
        if input_bam == NULL:
            sam_hdr_destroy(filtered_header)
            free(mapping.new_to_old_tid)
            free(mapping.old_to_new_tid)
            free(mapping)
            return -1
        opened_input_bam_here = True

    output_bam = hts_open(output_bam_path, b"wb")
    if output_bam == NULL:
        if opened_input_bam_here and input_bam != NULL:
            hts_close(input_bam)
        sam_hdr_destroy(filtered_header)
        free(mapping.new_to_old_tid)
        free(mapping.old_to_new_tid)
        free(mapping)
        return -1

    if num_threads > 1:
        hts_set_threads(output_bam, min_int32(num_threads, 4))

    hts_set_opt(input_bam, HTS_OPT_CACHE_SIZE, cache_size)
    hts_set_opt(output_bam, HTS_OPT_CACHE_SIZE, cache_size)

    if sam_hdr_write(output_bam, filtered_header) < 0:
        ret = -1
    else:
        clock_gettime(CLOCK_MONOTONIC, &ts_write_start)

        if existing_index != NULL:
            bam_index = existing_index
        else:
            bam_index = sam_index_load(input_bam, input_bam_path)
            if bam_index == NULL:
                ret = -1

        if ret == 0:
            record = bam_init1()
            if record == NULL:
                ret = -1
            else:
                for orig_tid in range(mapping.n_original_refs):
                    if ret != 0:
                        break

                    mapped_new = mapping.old_to_new_tid[orig_tid]
                    if mapped_new == <uint32_t>(-1):
                        continue

                    iterator = sam_itr_queryi(bam_index, orig_tid, 0, 0x7fffffff)
                    if iterator != NULL:
                        while sam_itr_next(input_bam, iterator, record) >= 0:
                            alignments_processed += 1

                            aln_read_length = record.core.l_qseq
                            if aln_read_length < min_read_length_c or aln_read_length > max_read_length_c:
                                continue

                            nm_tag_local = bam_aux_get(record, b"NM")
                            nm_val_local = -1
                            aln_ani = 0.0
                            if nm_tag_local != NULL:
                                nm_val_local = bam_aux2i(nm_tag_local)
                            if nm_val_local >= 0 and aln_read_length > 0:
                                aln_ani = (1.0 - (<double>nm_val_local / aln_read_length)) * 100.0
                            else:
                                aln_ani = 0.0

                            if aln_ani < min_read_ani_c:
                                continue

                            record.core.tid = <int32_t>mapped_new

                            if sam_write1(output_bam, filtered_header, record) < 0:
                                ret = -1
                                break

                            alignments_written += 1
                            ref_alignments += 1
                        hts_itr_destroy(iterator)
                        iterator = NULL

        if ret == 0:
            clock_gettime(CLOCK_MONOTONIC, &ts_write_end)
            write_sec = <double>(ts_write_end.tv_sec - ts_write_start.tv_sec) + <double>(ts_write_end.tv_nsec - ts_write_start.tv_nsec) / 1e9
            bf_nogil_logf_notime(NULL, "Filtered BAM complete: %lld/%lld alignments written to %s\n",
                                 alignments_written, alignments_processed, output_bam_path)
            bf_nogil_logf_notime(NULL, "Filtered BAM write time: %.6f sec\n", write_sec)

    if record != NULL:
        bam_destroy1(record)
    if existing_index == NULL and bam_index != NULL:
        hts_idx_destroy(bam_index)
    if opened_input_bam_here and input_bam != NULL:
        hts_close(input_bam)
    if output_bam != NULL:
        hts_close(output_bam)
    if filtered_header != NULL:
        sam_hdr_destroy(filtered_header)
    if mapping != NULL:
        if mapping.new_to_old_tid != NULL:
            free(mapping.new_to_old_tid)
        if mapping.old_to_new_tid != NULL:
            free(mapping.old_to_new_tid)
        free(mapping)

    return ret
