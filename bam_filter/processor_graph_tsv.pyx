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
"""
TSV writer extension for graph analysis output.

This module exposes a C-level function ``write_graph_tsv_c`` that other
Cython modules can cimport and call with low-level pointers. The function
writes a per-reference TSV (optionally gzipped) summarizing graph and
clustering metrics.
"""

from libc.stdlib cimport malloc, free, qsort
from libc.string cimport strlen
from libc.stdio cimport FILE, fopen, fclose, fprintf
from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t
from bam_filter.processor cimport MemoryPool
from bam_filter.processor_graph cimport ReferencePattern, ReferenceStats, MultPairExtended
from bam_filter.processor_mapping cimport ReferenceMapping

cdef extern from "htslib/sam.h":
    ctypedef struct sam_hdr_t
    const char* sam_hdr_tid2name(sam_hdr_t* header, int tid) nogil
    int64_t sam_hdr_tid2len(sam_hdr_t* header, int tid) nogil

cdef extern from "zlib.h":
    ctypedef void* gzFile
    gzFile gzopen(const char* path, const char* mode) nogil
    int gzclose(gzFile file) nogil
    int gzprintf(gzFile file, const char* format, ...) nogil


cdef inline bint c_str_endswith_gz(const char* s) nogil:
    """Check if C string ends with '.gz' extension."""
    if s == NULL:
        return 0
    cdef int L = <int>strlen(s)
    if L < 3:
        return 0
    if s[L-3] == 46 and s[L-2] == 103 and s[L-1] == 122:
        return 1
    return 0


cdef uint32_t UINT32_MAX = 0xFFFFFFFF

cdef int _multpair_extended_cmp_local(const void* a, const void* b) noexcept nogil:
    """Comparator for qsort: sort references by neighbor_count desc, then read_count desc."""
    cdef MultPairExtended* A = <MultPairExtended*>a
    cdef MultPairExtended* B = <MultPairExtended*>b

    if A.neighbor_count != B.neighbor_count:
        return -1 if A.neighbor_count > B.neighbor_count else 1

    if A.read_count != B.read_count:
        return -1 if A.read_count > B.read_count else 1

    return 0


cdef int write_graph_tsv_c(MemoryPool* pool, sam_hdr_t* bam_header,
                          ReferenceMapping* mapping, ReferencePattern* pattern_data,
                          ReferenceStats* ref_stats,
                          uint32_t* total_reads, uint32_t* multimap_reads, uint64_t* alignments_per_ref,
                          uint32_t* exact_connection_counts, double* co_mapping_averages,
                          uint64_t* max_co_mappings, uint32_t* co_mapping_counts,
                          double* neighbor_multimap_avg, double* neighbor_connections_avg,
                          uint32_t* neighbor_counts, uint32_t array_size,
                          double dataset_median_connections, int32_t min_read_count,
                          bint include_clustering, const char* tsv_path) noexcept nogil:
    """
    Write a TSV (optionally gzipped) summarizing graph analysis results per
    reference.

    This C-level function implements the TSV export previously located in
    :file:`bam_filter/processor_graph.pyx::write_graph_tsv` and was moved here
    to reduce the size of the main graph module. It is a low-level, ``nogil``
    function and must not perform Python API calls.

    Parameters
    ----------
    pool : MemoryPool*
        Pointer to the shared memory pool describing retained references and
        counts. Used for sizing and reference_count lookup.
    bam_header : sam_hdr_t*
        HTSlib header pointer used to map tid -> reference name/length. May
        be NULL, in which case reference name will be reported as ``unknown``.
    mapping : ReferenceMapping*
        Optional mapping object that translates new tids to original tids.
    pattern_data : ReferencePattern*
        Optional per-reference pattern/graph statistics produced by the
        clustering/graph pipeline. Used to populate leiden/community fields.
    ref_stats : ReferenceStats*
        Optional array of per-reference numeric statistics (read counts,
        score means, PMD aggregates). May be NULL when statistics are
        unavailable.
    total_reads, multimap_reads : uint32_t*
        Arrays of per-reference read counts. If provided, used to compute
        multimap percentages and per-reference read totals.
    alignments_per_ref : uint64_t*
        Array of total alignment counts per reference (optional).
    exact_connection_counts : uint32_t*
        If provided, exact neighbor counts for each reference (connection
        degrees) used for sorting and node-degree columns.
    co_mapping_averages : double*
        Per-reference average co-mappings per read (optional).
    max_co_mappings : uint64_t*
        Maximum co-mapping observed per reference (optional).
    co_mapping_counts : uint32_t*
        Optional counts used for distribution reporting.
    neighbor_multimap_avg, neighbor_connections_avg : double*
        Optional per-reference neighbor-level aggregates.
    neighbor_counts : uint32_t*
        Per-reference neighbor counts array (optional).
    array_size : uint32_t
        Length of the arrays passed above (safety bound for indexing).
    dataset_median_connections : double
        Median connections value for the dataset (used for contextual output).
    min_read_count : int32_t
        Minimum reads threshold used when deciding which references to emit.
    include_clustering : bint
        If true, Leiden/community columns will be filled when pattern_data is
        available.
    tsv_path : const char*
        Filesystem path for the TSV output. If it ends with ``.gz`` the file
        will be written gzipped using zlib's ``gzopen``/``gzprintf`` APIs.

    Returns
    -------
    int
        0 on success, -1 on any error (file open failure or allocation
        failure). This function performs internal allocations with ``malloc``
        and frees them before returning.

    Notes
    -----
    - This function is ``nogil`` and ``noexcept``. Do not introduce Python
      calls or raise exceptions inside it.
    - The caller is responsible for ensuring the arrays passed in are at
      least ``array_size`` long and that ``pool`` describes the correct
      reference_count.
    """

    cdef FILE* f = NULL
    cdef gzFile gzfp = NULL
    cdef uint32_t ref_idx, i
    cdef uint32_t node_degree
    cdef const char* ref_name_c
    cdef int32_t tid32
    cdef int64_t ref_len
    cdef double multimap_pct
    cdef uint32_t treads
    cdef uint32_t repeat_reads_val
    cdef uint32_t shared_reads_val
    cdef uint32_t connected_neighbors
    cdef double avg_co_mappings
    cdef uint64_t max_co_map
    cdef double neighbor_mm_rate
    cdef uint32_t unique_reads
    cdef uint32_t repeat_reads
    cdef uint32_t shared_reads
    cdef uint64_t align_count
    cdef double score_mean
    cdef double score_std
    cdef double score_min
    cdef double score_max
    cdef uint32_t pmd_available
    cdef double pmd_mean
    cdef double pmd_std
    cdef double pmd_min
    cdef double pmd_max
    cdef uint64_t pmd_nonzero_pct
    cdef uint32_t component_id, community_id
    cdef float community_cc
    cdef float individual_cc
    cdef float cc_threshold
    cdef char keep_flag
    cdef MultPairExtended* sorted_refs = NULL
    cdef uint32_t valid_count = 0

    cdef uint32_t UINT32_MAX = <uint32_t>0xFFFFFFFF

    if not tsv_path:
        return -1

    if c_str_endswith_gz(tsv_path):
        gzfp = gzopen(tsv_path, "wb")
        if gzfp == NULL:
            return -1

        gzprintf(gzfp, b"reference_name\treference_length_bp\tcomponent_id\tnode_degree\tleiden_community_id\tleiden_community_cc\tleiden_individual_cc\tleiden_cc_threshold\tleiden_keep_status\ttotal_reads\tunique_reads\trepeat_reads\tshared_reads\ttotal_alignments\tmultimap_percentage\tconnected_neighbors\tavg_co_mappings_per_read\tmax_co_mappings_observed\tneighbor_multimap_rate\talignment_score_mean\talignment_score_std\talignment_score_min\talignment_score_max\tpmd_available\tpmd_score_mean\tpmd_score_std\tpmd_score_min\tpmd_score_max\tpmd_nonzero_percentage\n")
    else:
        f = fopen(tsv_path, b"w")
        if not f:
            return -1

        fprintf(f, b"reference_name\treference_length_bp\tcomponent_id\tnode_degree\tleiden_community_id\tleiden_community_cc\tleiden_individual_cc\tleiden_cc_threshold\tleiden_keep_status\ttotal_reads\tunique_reads\trepeat_reads\tshared_reads\ttotal_alignments\tmultimap_percentage\tconnected_neighbors\tavg_co_mappings_per_read\tmax_co_mappings_observed\tneighbor_multimap_rate\talignment_score_mean\talignment_score_std\talignment_score_min\talignment_score_max\tpmd_available\tpmd_score_mean\tpmd_score_std\tpmd_score_min\tpmd_score_max\tpmd_nonzero_percentage\n")

    sorted_refs = <MultPairExtended*>malloc(pool.reference_count * sizeof(MultPairExtended))
    if not sorted_refs:
        if f: fclose(f)
        if gzfp != NULL: gzclose(gzfp)
        return -1

    for ref_idx in range(pool.reference_count):
        treads = 0
        if total_reads and ref_idx < pool.reference_count:
            treads = total_reads[ref_idx]
        if treads == 0:
            continue

        connected_neighbors = 0
        if ref_idx < array_size and exact_connection_counts:
            connected_neighbors = exact_connection_counts[ref_idx]

        sorted_refs[valid_count].neighbor_count = connected_neighbors
        sorted_refs[valid_count].read_count = treads
        sorted_refs[valid_count].idx = ref_idx

        if pattern_data and ref_idx < pool.reference_count:
            sorted_refs[valid_count].conc = <float>multimap_reads[ref_idx] / <float>total_reads[ref_idx] if total_reads[ref_idx] > 0 else 0.0
            sorted_refs[valid_count].net = pattern_data[ref_idx].graph.avg_comappings_per_read
        else:
            sorted_refs[valid_count].conc = 0.0
            sorted_refs[valid_count].net = 0.0

        valid_count += 1

    if valid_count > 1:
        qsort(sorted_refs, valid_count, sizeof(MultPairExtended), _multpair_extended_cmp_local)

    for i in range(valid_count):
        ref_idx = sorted_refs[i].idx

        ref_name_c = b"unknown"
        ref_len = 0

        if mapping and mapping.new_to_old_tid and ref_idx < mapping.n_retained_refs:
            tid32 = mapping.new_to_old_tid[ref_idx]
            if tid32 >= 0 and bam_header:
                ref_name_c = sam_hdr_tid2name(bam_header, tid32)
                ref_len = sam_hdr_tid2len(bam_header, tid32)

        if not ref_name_c:
            ref_name_c = b"unknown"

        treads = 0
        unique_reads = 0
        repeat_reads = 0
        shared_reads = 0
        if ref_stats and ref_idx < pool.reference_count:
            treads = ref_stats[ref_idx].total_reads
            unique_reads = ref_stats[ref_idx].unique_reads
            repeat_reads = ref_stats[ref_idx].repeat_reads
            shared_reads = ref_stats[ref_idx].shared_reads

        align_count = 0
        if ref_stats and ref_idx < pool.reference_count:
            align_count = ref_stats[ref_idx].alignment_count

        multimap_pct = 0.0
        if treads > 0:
            multimap_pct = 100.0 * <double>(repeat_reads + shared_reads) / <double>treads

        connected_neighbors = 0
        avg_co_mappings = 0.0
        max_co_map = 0
        if ref_idx < array_size:
            if exact_connection_counts:
                connected_neighbors = exact_connection_counts[ref_idx]
            if co_mapping_averages:
                avg_co_mappings = co_mapping_averages[ref_idx]
            elif pattern_data and ref_idx < pool.reference_count:
                avg_co_mappings = pattern_data[ref_idx].graph.avg_comappings_per_read

            if max_co_mappings:
                max_co_map = max_co_mappings[ref_idx]
            elif pattern_data and ref_idx < pool.reference_count:
                max_co_map = pattern_data[ref_idx].graph.max_comappings

        neighbor_mm_rate = 0.0
        if ref_idx < array_size and neighbor_multimap_avg:
            neighbor_mm_rate = neighbor_multimap_avg[ref_idx] * 100.0

        score_mean = 0.0
        score_std = 0.0
        score_min = 0.0
        score_max = 0.0
        if ref_stats and ref_idx < pool.reference_count:
            score_mean = ref_stats[ref_idx].score_mean
            score_std = ref_stats[ref_idx].score_std
            score_min = ref_stats[ref_idx].score_min
            score_max = ref_stats[ref_idx].score_max

        pmd_available = 0
        pmd_mean = 0.0
        pmd_std = 0.0
        pmd_min = 0.0
        pmd_max = 0.0
        pmd_nonzero_pct = 0
        if ref_stats and ref_idx < pool.reference_count:
            pmd_available = 1 if ref_stats[ref_idx].pmd_available else 0
            pmd_mean = ref_stats[ref_idx].pmd_mean
            pmd_std = ref_stats[ref_idx].pmd_std
            pmd_min = ref_stats[ref_idx].pmd_min
            pmd_max = ref_stats[ref_idx].pmd_max
            if align_count > 0:
                pmd_nonzero_pct = <uint64_t>(100.0 * <double>ref_stats[ref_idx].pmd_nonzero_count / <double>align_count)

        if avg_co_mappings != avg_co_mappings or avg_co_mappings > 1e9 or avg_co_mappings < -1e9:
            avg_co_mappings = 0.0

        component_id = <uint32_t>UINT32_MAX
        community_id = <uint32_t>UINT32_MAX
        community_cc = 0.0
        individual_cc = 0.0
        cc_threshold = 0.0
        keep_flag = 0

        node_degree = 0
        if pattern_data and ref_idx < pool.reference_count:
            component_id = pattern_data[ref_idx].component_id
            node_degree = pattern_data[ref_idx].node_degree
            community_id = pattern_data[ref_idx].leiden_community_id
            community_cc = pattern_data[ref_idx].leiden_community_cc
            individual_cc = pattern_data[ref_idx].leiden_individual_cc
            cc_threshold = pattern_data[ref_idx].leiden_cc_threshold
            keep_flag = pattern_data[ref_idx].leiden_keep_flag

        if gzfp != NULL:
            if community_id == UINT32_MAX:
                if keep_flag:
                    gzprintf(gzfp, b"%s\t%ld\tcomp_%u\t%u\tsingleton\t%.6f\t%.6f\t%.6f\tkept\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\n",
                        ref_name_c, ref_len, component_id, node_degree, community_cc, individual_cc, cc_threshold,
                        treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                        connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                        score_mean, score_std, score_min, score_max, pmd_available,
                        pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct)
                else:
                    gzprintf(gzfp, b"%s\t%ld\tcomp_%u\t%u\tsingleton\t%.6f\t%.6f\t%.6f\tremoved\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\n",
                        ref_name_c, ref_len, component_id, node_degree, community_cc, individual_cc, cc_threshold,
                        treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                        connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                        score_mean, score_std, score_min, score_max, pmd_available,
                        pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct)
            else:
                if keep_flag:
                    gzprintf(gzfp, b"%s\t%ld\tcomp_%u\t%u\tcomm_%u\t%.6f\t%.6f\t%.6f\tkept\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\n",
                        ref_name_c, ref_len, component_id, node_degree, community_id, community_cc, individual_cc, cc_threshold,
                        treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                        connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                        score_mean, score_std, score_min, score_max, pmd_available,
                        pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct)
                else:
                    gzprintf(gzfp, b"%s\t%ld\tcomp_%u\t%u\tcomm_%u\t%.6f\t%.6f\t%.6f\tremoved\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\n",
                        ref_name_c, ref_len, component_id, node_degree, community_id, community_cc, individual_cc, cc_threshold,
                        treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                        connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                        score_mean, score_std, score_min, score_max, pmd_available,
                        pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct)
        else:
            if community_id == UINT32_MAX:
                if keep_flag:
                    fprintf(f, b"%s\t%ld\tcomp_%u\t%u\tsingleton\t%.6f\t%.6f\t%.6f\tkept\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\n",
                        ref_name_c, ref_len, component_id, node_degree, community_cc, individual_cc, cc_threshold,
                        treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                        connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                        score_mean, score_std, score_min, score_max, pmd_available,
                        pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct)
                else:
                    fprintf(f, b"%s\t%ld\tcomp_%u\t%u\tsingleton\t%.6f\t%.6f\t%.6f\tremoved\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\n",
                        ref_name_c, ref_len, component_id, node_degree, community_cc, individual_cc, cc_threshold,
                        treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                        connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                        score_mean, score_std, score_min, score_max, pmd_available,
                        pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct)
            else:
                if keep_flag:
                    fprintf(f, b"%s\t%ld\tcomp_%u\t%u\tcomm_%u\t%.6f\t%.6f\t%.6f\tkept\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\n",
                        ref_name_c, ref_len, component_id, node_degree, community_id, community_cc, individual_cc, cc_threshold,
                        treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                        connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                        score_mean, score_std, score_min, score_max, pmd_available,
                        pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct)
                else:
                    fprintf(f, b"%s\t%ld\tcomp_%u\t%u\tcomm_%u\t%.6f\t%.6f\t%.6f\tremoved\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\n",
                        ref_name_c, ref_len, component_id, node_degree, community_id, community_cc, individual_cc, cc_threshold,
                        treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                        connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                        score_mean, score_std, score_min, score_max, pmd_available,
                        pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct)

    free(sorted_refs)
    if gzfp != NULL:
        gzclose(gzfp)
    if f != NULL:
        fclose(f)
    return 0
