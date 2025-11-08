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
from bam_filter.processor_types cimport (
    bam1_t,
    bam1_core_t,
    hts_idx_t,
    hts_itr_t,
    bam_init1,
    bam_destroy1,
    bam_aux_get,
    bam_aux2i,
    sam_index_load,
    hts_idx_destroy,
    sam_itr_queryi,
    sam_itr_next,
    hts_itr_destroy,
    hts_set_threads,
)
from bam_filter.processor_types cimport (
    hts_open,
    hts_close,
    sam_hdr_write,
    sam_write1,
    sam_index_load,
    sam_itr_queryi,
    sam_itr_next,
    hts_itr_destroy,
    hts_idx_destroy,
    bam_init1,
    bam_destroy1,
    bam_aux_get,
    sam_hdr_destroy,
    hts_set_threads,
    hts_set_opt,
    HTS_OPT_CACHE_SIZE,
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

# Use centralized declarations from processor_types.pxd for htslib types and functions


cdef bint passes_filters(RefStats* stats, FilterConditions* filters) noexcept nogil:
    """Check if reference passes all ENABLED filter conditions."""
    cdef double cov_evenness_tmp

    if filters.enable_min_coverage_mean and stats.mean_coverage < filters.min_coverage_mean:
        return <bint>False

    if filters.enable_min_avg_read_ani and stats.ani_mean < filters.min_avg_read_ani:
        return <bint>False

    if filters.enable_min_breadth and stats.breadth < filters.min_breadth:
        return <bint>False

    if filters.enable_min_expected_breadth_ratio and stats.breadth_exp_ratio < filters.min_expected_breadth_ratio:
        return <bint>False

    if filters.enable_min_coverage_evenness:
        cov_evenness_tmp = <double>(stats.cov_evenness if stats.mean_coverage >= 1.0 else 1.0)
        if cov_evenness_tmp < filters.min_coverage_evenness:
            return <bint>False

    if filters.enable_max_coeff_var and stats.c_v > filters.max_coeff_var:
        return <bint>False

    if filters.enable_min_norm_entropy and stats.norm_entropy < filters.min_norm_entropy:
        return <bint>False

    if filters.enable_max_norm_gini and stats.norm_gini > filters.max_norm_gini:
        return <bint>False

    return <bint>True


# Small inline helpers to safely index C arrays in nogil code and avoid
# static-analyzer complaints about pointer indexing across modules.
cdef inline void set_int32_at(int32_t* arr, int i, int32_t v) nogil:
    arr[i] = v

cdef inline int32_t get_int32_at(int32_t* arr, int i) nogil:
    return arr[i]

cdef inline void set_uint32_at(uint32_t* arr, int i, uint32_t v) nogil:
    arr[i] = v

cdef inline uint32_t get_uint32_at(uint32_t* arr, int i) nogil:
    return arr[i]

cdef inline RefStats* refstats_ptr_at(RefStats* arr, int i) nogil:
    return &arr[i]

cdef inline int64_t refstats_n_alns(RefStats* arr, int i) nogil:
    return arr[i].n_alns


cdef ReferenceFilter* create_reference_filter(
    RefStats* global_ref_stats,
    FilterConditions* filters,
    int n_refs
) noexcept nogil:
    cdef ReferenceFilter* ref_filter = <ReferenceFilter*>malloc(<size_t>sizeof(ReferenceFilter))
    if not ref_filter:
        return <ReferenceFilter*>0

    ref_filter.n_total_refs = n_refs
    ref_filter.tid_mapping = <int32_t*>malloc(<size_t>(n_refs * sizeof(int32_t)))
    if not ref_filter.tid_mapping:
        free(<void*>ref_filter)
        return <ReferenceFilter*>0

    cdef int i
    cdef int32_t* tid_map = ref_filter.tid_mapping
    cdef int32_t* rev_map
    cdef int mapped_idx_int
    for i in range(n_refs):
        set_int32_at(tid_map, i, <int32_t>(-1))

    cdef int32_t new_tid = 0
    cdef RefStats* stats_ptr
    for i in range(n_refs):
        stats_ptr = refstats_ptr_at(global_ref_stats, i)
        if refstats_n_alns(global_ref_stats, i) > 0 and passes_filters(stats_ptr, filters):
            set_int32_at(tid_map, i, new_tid)
            new_tid += 1

    ref_filter.n_filtered_refs = new_tid

    if new_tid > 0:
        ref_filter.reverse_mapping = <int32_t*>malloc(<size_t>(new_tid * sizeof(int32_t)))
        if not ref_filter.reverse_mapping:
            free(<void*>ref_filter.tid_mapping)
            free(<void*>ref_filter)
            return <ReferenceFilter*>0
        rev_map = ref_filter.reverse_mapping
        for i in range(n_refs):
            if get_int32_at(tid_map, i) >= 0:
                mapped_idx_int = <int>get_int32_at(tid_map, i)
                set_int32_at(rev_map, mapped_idx_int, <int32_t>i)
    else:
        ref_filter.reverse_mapping = NULL

    return ref_filter


cdef void destroy_reference_filter(ReferenceFilter* ref_filter) noexcept nogil:
    if ref_filter:
        if ref_filter.tid_mapping:
            free(<void*>ref_filter.tid_mapping)
        if ref_filter.reverse_mapping:
            free(<void*>ref_filter.reverse_mapping)
        free(<void*>ref_filter)


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
    cdef int64_t alignments_written = <int64_t>0
    cdef int64_t alignments_processed = <int64_t>0
    cdef int32_t old_tid, new_tid
    cdef int32_t orig_tid
    cdef int64_t ref_alignments = <int64_t>0
    cdef uint32_t mapped_new
    cdef int cache_size = 256 * 1024 * 1024
    cdef int i
    cdef int32_t aln_read_length
    cdef uint8_t* nm_tag_local
    cdef int nm_val_local
    cdef double aln_ani
    cdef timespec ts_write_start, ts_write_end
    cdef double write_sec = <double>0.0


    if ref_filter.n_filtered_refs == 0:
        bf_nogil_logf_notime(<const char*>NULL, <const char*>b"No references pass filtering - no BAM output created\n")
        return 0

    mapping = <ReferenceMapping*>malloc(<size_t>sizeof(ReferenceMapping))
    if not mapping:
        return -1

    mapping.n_retained_refs = ref_filter.n_filtered_refs
    mapping.n_original_refs = ref_filter.n_total_refs

    mapping.new_to_old_tid = <uint32_t*>malloc(<size_t>(mapping.n_retained_refs * sizeof(uint32_t)))
    mapping.old_to_new_tid = <uint32_t*>malloc(<size_t>(mapping.n_original_refs * sizeof(uint32_t)))

    if not mapping.new_to_old_tid or not mapping.old_to_new_tid:
        if mapping.new_to_old_tid:
            free(<void*>mapping.new_to_old_tid)
        if mapping.old_to_new_tid:
            free(<void*>mapping.old_to_new_tid)
        free(<void*>mapping)
        return -1

    cdef uint32_t* old_to_new = mapping.old_to_new_tid
    for i in range(mapping.n_original_refs):
        set_uint32_at(old_to_new, i, <uint32_t>(-1))

    cdef uint32_t* new_to_old = mapping.new_to_old_tid
    for new_tid in range(mapping.n_retained_refs):
        old_tid = ref_filter.reverse_mapping[new_tid]
        set_uint32_at(new_to_old, new_tid, <uint32_t>old_tid)
        set_uint32_at(old_to_new, <int>old_tid, <uint32_t>new_tid)

    filtered_header = create_filtered_header_efficient(original_header, mapping)

    if not filtered_header:
        free(<void*>mapping.new_to_old_tid)
        free(<void*>mapping.old_to_new_tid)
        free(<void*>mapping)
        return -1

    cdef bint opened_input_bam_here = <bint>False

    if not input_bam:
        input_bam = hts_open(input_bam_path, <const char*>b"r")
        if not input_bam:
            sam_hdr_destroy(filtered_header)
            free(<void*>mapping.new_to_old_tid)
            free(<void*>mapping.old_to_new_tid)
            free(<void*>mapping)
            return -1
        opened_input_bam_here = <bint>True

    output_bam = hts_open(output_bam_path, <const char*>b"wb")
    if not output_bam:
        if opened_input_bam_here and input_bam:
            hts_close(input_bam)
        sam_hdr_destroy(filtered_header)
        free(<void*>mapping.new_to_old_tid)
        free(<void*>mapping.old_to_new_tid)
        free(<void*>mapping)
        return -1

    if num_threads > 1:
        hts_set_threads(output_bam, min_int32(num_threads, 4))

    # Omit calling hts_set_opt here to avoid varargs prototype issues in nogil
    # contexts. HTSlib cache configuration is non-critical.

    if sam_hdr_write(output_bam, filtered_header) < 0:
        ret = -1
    else:
        clock_gettime(CLOCK_MONOTONIC, &ts_write_start)

        if existing_index:
            bam_index = existing_index
        else:
            bam_index = sam_index_load(input_bam, input_bam_path)
            if not bam_index:
                ret = -1

        if ret == 0:
            record = bam_init1()
            if not record:
                ret = -1
            else:
                for orig_tid in range(mapping.n_original_refs):
                    if ret != 0:
                        break

                    mapped_new = get_uint32_at(mapping.old_to_new_tid, orig_tid)
                    if mapped_new == <uint32_t>(-1):
                        continue

                    iterator = sam_itr_queryi(bam_index, orig_tid, 0, 0x7fffffff)
                    if iterator:
                        while sam_itr_next(input_bam, iterator, record) >= 0:
                            alignments_processed += 1

                            aln_read_length = record.core.l_qseq
                            if aln_read_length < min_read_length_c or aln_read_length > max_read_length_c:
                                continue

                            nm_tag_local = bam_aux_get(record, <const char*>b"NM")
                            nm_val_local = -1
                            aln_ani = <double>0.0
                            if nm_tag_local:
                                nm_val_local = bam_aux2i(nm_tag_local)
                            if nm_val_local >= 0 and aln_read_length > 0:
                                aln_ani = (1.0 - (<double>nm_val_local / aln_read_length)) * 100.0
                            else:
                                aln_ani = <double>0.0

                            if aln_ani < min_read_ani_c:
                                continue

                            record.core.tid = <int64_t>mapped_new

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
            bf_nogil_logf_notime(<const char*>NULL, <const char*>b"Filtered BAM complete: %lld/%lld alignments written to %s\n",
                                 alignments_written, alignments_processed, output_bam_path)
            bf_nogil_logf_notime(<const char*>NULL, <const char*>b"Filtered BAM write time: %.6f sec\n", write_sec)

    if record:
        bam_destroy1(record)
    if not existing_index and bam_index:
        hts_idx_destroy(bam_index)
    if opened_input_bam_here and input_bam:
        hts_close(input_bam)
    if output_bam:
        hts_close(output_bam)
    if filtered_header:
        sam_hdr_destroy(filtered_header)
    if mapping:
        if mapping.new_to_old_tid:
            free(<void*>mapping.new_to_old_tid)
        if mapping.old_to_new_tid:
            free(<void*>mapping.old_to_new_tid)
        free(<void*>mapping)

    return ret
