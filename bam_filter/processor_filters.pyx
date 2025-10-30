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

"""Alignment filtering strategies.

Implements probability-based and cluster-aware filtering for removing
low-confidence alignments and spurious references based on graph topology.
"""

from cython.parallel cimport prange, threadid

from libc.stdlib cimport calloc, free, malloc, realloc
from libc.string cimport memset
from libc.stdint cimport int64_t, int32_t, uint32_t, uint64_t
from libc.math cimport log2, fmin, fmax, exp, sqrt, pow, INFINITY

cdef uint32_t UINT32_MAX = 0xFFFFFFFF
from bam_filter.processor cimport MemoryPool
from bam_filter.processor_graph cimport ReferencePattern, ReferenceStats, ReadIndex
from bam_filter.processor_types cimport EMAlgorithmConfig, PrecomputedWeights
from bam_filter.processor_precomputed cimport (
    create_precomputed_weights, 
    free_precomputed_weights, update_precomputed_weights
)
from bam_filter.processor_em cimport get_reference_weights
from bam_filter.processor_fast_math cimport stable_log_sum_exp

from bam_filter.processor_graph_ops cimport (
    WeightedGraph, GraphNode,
    create_weighted_graph, destroy_weighted_graph, add_edge,
    build_weighted_graph_from_alignments, prune_low_weight_edges,
    calculate_graph_statistics
)
from bam_filter.processor_leiden_igraph cimport (
    leiden_clustering,
    LeidenResults,
)
from bam_filter.processor_igraph cimport *
from bam_filter.processor_graph cimport write_graph_tsv
from bam_filter.processor_graph_export cimport export_graph_graphml
from bam_filter.processor_mapping cimport ReferenceMapping
cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil
cdef extern from "htslib/sam.h":
    ctypedef struct sam_hdr_t


cdef int apply_probability_filtering(MemoryPool* pool, EMAlgorithmConfig* config) except -1 nogil:
    """Filter alignments by EM-derived posterior probability.

    Removes alignments below minimum probability threshold and/or fraction of
    maximum probability per read. Uses precomputed reference weights to calculate
    posteriors efficiently. Compacts alignment array and ZP values in-place.

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool containing alignments
    config : EMAlgorithmConfig*
        Configuration with probability thresholds and thread count

    Returns
    -------
    int
        0 on success, -1 on error

    Notes
    -----
    PMD data is stored in alignment structures and moves automatically during
    compaction. Rebuilds read indexing structures after filtering.
    """
    cdef uint32_t read_idx, ref_idx, rid
    cdef int64_t start_pos, end_pos, ai
    cdef uint32_t alignment_count
    cdef float alignment_score, uniform_zp
    cdef double log_lik 
    cdef double log_weighted, log_norm, posterior
    cdef double NEG_INF = -1e20
    cdef int64_t alignments_removed
    cdef int64_t write_idx = 0
    cdef bint keep_alignment
    cdef double min_threshold = fmax(config.minimum_probability_threshold, 1e-8)
    cdef double fraction_threshold = config.probability_fraction_filter
    cdef double* reference_weights = get_reference_weights(pool)
    cdef PrecomputedWeights* precomp = NULL
    cdef float* read_max_probs
    cdef int32_t* survivors_per_read
    cdef int64_t survivors_total = 0
    cdef uint64_t start_pos_rebuild
    cdef double thr
    cdef double p

    precomp = create_precomputed_weights(pool.reference_count)
    if not precomp:
        return -1
    update_precomputed_weights(precomp, reference_weights)

    bf_nogil_logf_notime(
        b"FILTER",
        "probability_filter: start pmd_output=%s total_alignments=%lld",
        b"enabled" if pool.pmd_enabled_for_output else b"disabled",
        <long long>pool.alignment_count,
    )
    if pool.scratch_read_max_probs == NULL or pool.scratch_survivors_per_read == NULL or pool.scratch_unique_read_count < <int32_t>pool.unique_read_count:
        if pool.scratch_read_max_probs != NULL:
            free(pool.scratch_read_max_probs)
        if pool.scratch_survivors_per_read != NULL:
            free(pool.scratch_survivors_per_read)
        pool.scratch_read_max_probs = <float*>calloc(pool.unique_read_count, sizeof(float))
        pool.scratch_survivors_per_read = <int32_t*>calloc(pool.unique_read_count, sizeof(int32_t))
        pool.scratch_unique_read_count = pool.unique_read_count
    else:
        memset(pool.scratch_read_max_probs, 0, pool.unique_read_count * sizeof(float))
        memset(pool.scratch_survivors_per_read, 0, pool.unique_read_count * sizeof(int32_t))
    read_max_probs = pool.scratch_read_max_probs
    survivors_per_read = pool.scratch_survivors_per_read
    for read_idx in prange(pool.unique_read_count, nogil=True, schedule='static', num_threads=config.thread_count):
        alignment_count = pool.read_alignment_counts[read_idx]
        if alignment_count == 0:
            continue
        start_pos = pool.read_alignment_starts[read_idx]
        end_pos = start_pos + alignment_count
        if alignment_count == 1:
            read_max_probs[read_idx] = 1.0
            if 1.0 >= min_threshold and (fraction_threshold == 0.0 or 1.0 >= fraction_threshold * 1.0):
                survivors_per_read[read_idx] = 1
            continue

        log_norm = NEG_INF
        for ai in range(start_pos, end_pos):
            ref_idx = pool.alignments[ai].reference_index
            if ref_idx < pool.reference_count:
                alignment_score = pool.alignments[ai].alignment_score
                log_lik = <double>alignment_score
                log_weighted = precomp.log_weights[ref_idx] + log_lik
                log_norm = stable_log_sum_exp(log_norm, log_weighted)
        if log_norm == NEG_INF:
            uniform_zp = 1.0 / alignment_count if alignment_count > 0 else 0.0
            read_max_probs[read_idx] = uniform_zp
            if uniform_zp >= min_threshold and (fraction_threshold == 0.0 or uniform_zp >= fraction_threshold * uniform_zp):
                survivors_per_read[read_idx] = alignment_count
            continue

        p = 0.0
        for ai in range(start_pos, end_pos):
            ref_idx = pool.alignments[ai].reference_index
            if ref_idx < pool.reference_count:
                alignment_score = pool.alignments[ai].alignment_score
                log_lik = <double>alignment_score
                log_weighted = precomp.log_weights[ref_idx] + log_lik
                posterior = exp(log_weighted - log_norm)
                if posterior > p:
                    p = posterior
        read_max_probs[read_idx] = <float>p
        thr = fraction_threshold * p
        for ai in range(start_pos, end_pos):
            ref_idx = pool.alignments[ai].reference_index
            if ref_idx < pool.reference_count:
                alignment_score = pool.alignments[ai].alignment_score
                log_lik = <double>alignment_score
                log_weighted = precomp.log_weights[ref_idx] + log_lik
                posterior = exp(log_weighted - log_norm)
                if posterior >= min_threshold and (fraction_threshold == 0.0 or posterior >= thr):
                    survivors_per_read[read_idx] += 1

    survivors_total = 0
    for read_idx in range(pool.unique_read_count):
        survivors_total += survivors_per_read[read_idx]

    bf_nogil_logf_notime(
        b"FILTER",
        "probability_filter: survivors=%lld total=%lld",
        <long long>survivors_total,
        <long long>pool.alignment_count,
    )

    if survivors_total <= 0:
        bf_nogil_logf_notime(b"WARN", "probability_filter: survivors=0 (aborting)")
        free_precomputed_weights(precomp)
        return -1
    if pool.precomputed_zp_values:
        free(pool.precomputed_zp_values)
    pool.precomputed_zp_values = <float*>malloc(survivors_total * sizeof(float))
    if not pool.precomputed_zp_values:
        free_precomputed_weights(precomp)
        return -1

    write_idx = 0
    for read_idx in range(pool.unique_read_count):
        alignment_count = pool.read_alignment_counts[read_idx]
        if alignment_count == 0:
            continue

        start_pos = pool.read_alignment_starts[read_idx]
        end_pos = start_pos + alignment_count

        if alignment_count == 1:
            p = 1.0
            if p >= min_threshold and (fraction_threshold == 0.0 or p >= fraction_threshold * 1.0):
                if write_idx != start_pos:
                    pool.alignments[write_idx] = pool.alignments[start_pos]
                pool.precomputed_zp_values[write_idx] = 1.0
                write_idx += 1
            continue

        log_norm = NEG_INF
        for ai in range(start_pos, end_pos):
            ref_idx = pool.alignments[ai].reference_index
            if ref_idx < pool.reference_count:
                alignment_score = pool.alignments[ai].alignment_score
                log_weighted = precomp.log_weights[ref_idx] + alignment_score
                log_norm = stable_log_sum_exp(log_norm, log_weighted)

        if log_norm == NEG_INF:
            uniform_zp = 1.0 / alignment_count if alignment_count > 0 else 0.0
            if uniform_zp >= min_threshold and (fraction_threshold == 0.0 or uniform_zp >= fraction_threshold * uniform_zp):
                for ai in range(start_pos, end_pos):
                    if write_idx != ai:
                        pool.alignments[write_idx] = pool.alignments[ai]
                    pool.precomputed_zp_values[write_idx] = uniform_zp
                    write_idx += 1
            continue

        thr = (<double>fraction_threshold) * (<double>read_max_probs[read_idx])
        for ai in range(start_pos, end_pos):
            ref_idx = pool.alignments[ai].reference_index
            if ref_idx < pool.reference_count:
                alignment_score = pool.alignments[ai].alignment_score
                log_weighted = precomp.log_weights[ref_idx] + alignment_score
                posterior = exp(log_weighted - log_norm)
                keep_alignment = (posterior >= min_threshold) and (fraction_threshold == 0.0 or posterior >= thr)
                if keep_alignment:
                    if write_idx != ai:
                        pool.alignments[write_idx] = pool.alignments[ai]
                    if posterior < 1e-12:
                        posterior = 1e-12
                    elif posterior > 0.999:
                        posterior = 0.999
                    pool.precomputed_zp_values[write_idx] = <float>posterior
                    write_idx += 1

    alignments_removed = pool.alignment_count - write_idx
    pool.alignment_count = write_idx
    pool.zp_values_computed = True

    bf_nogil_logf_notime(
        b"FILTER",
        "probability_filter: compacted_alignments=%llu removed=%lld",
        <unsigned long long>write_idx,
        <long long>alignments_removed,
    )

    memset(pool.read_alignment_counts, 0, pool.unique_read_count * sizeof(uint32_t))
    for ai in range(pool.alignment_count):
        rid = pool.alignments[ai].read_index
        if rid < pool.unique_read_count:
            pool.read_alignment_counts[rid] += 1
        else:
            bf_nogil_logf_notime(
                NULL,
                "probability_filter_error: invalid_read_id=%u alignment=%llu",
                rid,
                <unsigned long long>ai,
            )

    start_pos_rebuild = 0
    for rid in range(pool.unique_read_count):
        pool.read_alignment_starts[rid] = start_pos_rebuild
        start_pos_rebuild += pool.read_alignment_counts[rid]

    pool.final_unique_reads = 0
    for rid in range(pool.unique_read_count):
        if pool.read_alignment_counts[rid] > 0:
            pool.final_unique_reads += 1

    free_precomputed_weights(precomp)
    return 0


cdef int apply_reference_filtering(MemoryPool* pool, int32_t min_read_count) noexcept nogil:
    """Filter references below minimum read count threshold.

    Removes all alignments to references with fewer than min_read_count alignments.
    Compacts alignment array and ZP values, rebuilds read indexing structures.

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool containing alignments
    min_read_count : int32_t
        Minimum read count threshold

    Returns
    -------
    int
        0 on success, -1 on error

    Notes
    -----
    PMD data is integrated in alignment structures and moves automatically.
    """
    cdef int64_t* reference_counts = <int64_t*>calloc(pool.reference_count, sizeof(int64_t))
    cdef char* reference_keep_flag = <char*>calloc(pool.reference_count, sizeof(char))
    cdef int64_t alignment_idx, new_alignment_idx = 0
    cdef uint32_t ref_id
    cdef int64_t valid_references = 0
    cdef int64_t original_count = pool.alignment_count
    cdef uint32_t current_rid
    cdef int64_t count_errors = 0
    cdef uint64_t start_pos = 0
    cdef uint32_t rid
    cdef float* new_zp_values = NULL

    if not reference_counts or not reference_keep_flag:
        if reference_counts: free(reference_counts)
        if reference_keep_flag: free(reference_keep_flag)
        return -1

    for alignment_idx in range(pool.alignment_count):
        ref_id = pool.alignments[alignment_idx].reference_index
        if ref_id < pool.reference_count:
            reference_counts[ref_id] += 1
    cdef int64_t zero_count = 0
    cdef int64_t below_min_count = 0
    for ref_id in range(pool.reference_count):
        if reference_counts[ref_id] >= min_read_count:
            reference_keep_flag[ref_id] = 1
            valid_references += 1
        else:
            reference_keep_flag[ref_id] = 0
            if reference_counts[ref_id] == 0:
                zero_count += 1
            else:
                below_min_count += 1

    bf_nogil_logf_notime(
        b"FILTER-PMD",
        "%ld/%u references pass min_read_count=%d\n",
        valid_references,
        pool.reference_count,
        min_read_count,
    )
    bf_nogil_logf_notime(
        b"FILTER-PMD",
        "  Diagnostic: zero-count refs=%ld, refs_with_reads_below_min=%ld, min_read_count=%d\n",
        zero_count,
        below_min_count,
        min_read_count,
    )
    cdef int64_t surviving_alignments = 0
    for alignment_idx in range(pool.alignment_count):
        ref_id = pool.alignments[alignment_idx].reference_index
        if ref_id < pool.reference_count and reference_keep_flag[ref_id]:
            surviving_alignments += 1

    if surviving_alignments > 0:
        new_zp_values = <float*>malloc(surviving_alignments * sizeof(float))
        if not new_zp_values:
            free(reference_counts)
            free(reference_keep_flag)
            return -1
    else:
        bf_nogil_logf_notime(b"FILTER-PMD", "ERROR: No alignments survive reference filtering\n")
        free(reference_counts)
        free(reference_keep_flag)
        return -1

    new_alignment_idx = 0
    for alignment_idx in range(pool.alignment_count):
        ref_id = pool.alignments[alignment_idx].reference_index

        if (ref_id < pool.reference_count and reference_keep_flag[ref_id]):
            if new_alignment_idx != alignment_idx:
                pool.alignments[new_alignment_idx] = pool.alignments[alignment_idx]

            if new_zp_values and pool.precomputed_zp_values:
                new_zp_values[new_alignment_idx] = pool.precomputed_zp_values[alignment_idx]

            new_alignment_idx += 1
    if pool.precomputed_zp_values:
        free(pool.precomputed_zp_values)
    pool.precomputed_zp_values = new_zp_values

    cdef int64_t alignments_removed = pool.alignment_count - new_alignment_idx
    pool.alignment_count = new_alignment_idx

    bf_nogil_logf_notime(
        b"FILTER-PMD",
        "Removed %ld alignments, %ld remain\n",
        alignments_removed,
        pool.alignment_count,
    )
    bf_nogil_logf_notime(b"FILTER-PMD", "  All data (including integrated PMD) moved together\n")

    if alignments_removed > 0:
        bf_nogil_logf_notime(b"FILTER-PMD", "Rebuilding read indices...\n")

        memset(pool.read_alignment_counts, 0, pool.unique_read_count * sizeof(uint32_t))
        for alignment_idx in range(pool.alignment_count):
            current_rid = pool.alignments[alignment_idx].read_index
            if current_rid < pool.unique_read_count:
                pool.read_alignment_counts[current_rid] += 1
            else:
                count_errors += 1

        if count_errors > 0:
            bf_nogil_logf_notime(b"FILTER-PMD", "ERROR: Found %ld invalid read indices during rebuild\n", count_errors)
            free(reference_counts)
            free(reference_keep_flag)
            return -1

        for rid in range(pool.unique_read_count):
            pool.read_alignment_starts[rid] = start_pos
            start_pos += pool.read_alignment_counts[rid]

        if start_pos != <uint64_t>pool.alignment_count:
            bf_nogil_logf_notime(
                b"FILTER-PMD",
                "ERROR: Rebuild mismatch: calculated %lu, expected %ld\n",
                <unsigned long>start_pos,
                pool.alignment_count,
            )
            free(reference_counts)
            free(reference_keep_flag)
            return -1

        bf_nogil_logf_notime(b"FILTER-PMD", "Read index rebuild completed\n")

    free(reference_counts)
    free(reference_keep_flag)

    bf_nogil_logf_notime(b"FILTER-PMD", "Completed - all data synchronized\n")
    return 0


cdef int apply_cluster_aware_filtering(MemoryPool* pool,
                                       ReferencePattern* pattern_data,
                                       ReferenceStats* ref_stats,
                                       ReadIndex* read_index,
                                       uint32_t array_size,
                                       int32_t min_read_count,
                                       float score_threshold,
                                       bint verbose,
                                       int use_leiden,
                                       double leiden_resolution,
                                       bint leiden_parallel,
                                       int leiden_max_iterations,
                                       uint32_t graph_min_edge_weight,
                                       int32_t thread_count,
                                       int32_t iforest_n_trees,
                                       uint32_t iforest_subsample_size,
                                       double iforest_contamination,
                                       uint32_t iforest_random_seed,
                                       uint32_t lof_k,
                                       double lof_contamination,
                                       double zscore_threshold,
                                       WeightedGraph* existing_graph,
                                       sam_hdr_t* bam_header,
                                       ReferenceMapping* mapping,
                                       const char* tsv_export_path,
                                       const char* graph_export_path,
                                       int outlier_method) except -1 nogil:
    """Apply cluster-aware filtering based on graph topology and Leiden clustering.

    Uses igraph for graph construction, component detection, Leiden clustering,
    and weighted clustering coefficient calculation. Applies broken stick model
    to identify and remove promiscuous hub references within communities.

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool containing alignments
    pattern_data : ReferencePattern*
        Array for storing clustering results (output)
    ref_stats : ReferenceStats*
        Reference statistics array
    read_index : ReadIndex*
        Read indexing structure (unused, retained for API compatibility)
    array_size : uint32_t
        Size of pattern_data and ref_stats arrays
    min_read_count : int32_t
        Minimum read count filter threshold
    score_threshold : float
        Score threshold (unused, retained for API compatibility)
    verbose : bint
        Enable verbose logging
    use_leiden : int
        Enable Leiden clustering (currently always enabled)
    leiden_resolution : double
        Leiden resolution parameter
    leiden_parallel : bint
        Enable parallel Leiden (currently unused)
    leiden_max_iterations : int
        Maximum Leiden iterations
    graph_min_edge_weight : uint32_t
        Minimum edge weight for graph construction
    thread_count : int32_t
        Thread count for parallel processing
    existing_graph : WeightedGraph*
        Pre-built graph (required)
    bam_header : sam_hdr_t*
        BAM header for TSV/GraphML export
    mapping : ReferenceMapping*
        Reference mapping for export
    tsv_export_path : const char*
        Path for TSV output (NULL to skip)
    graph_export_path : const char*
        Path for GraphML output (NULL to skip)

    Returns
    -------
    int
        0 on success, -1 on error

    Notes
    -----
    Workflow:
    1. Verify pre-built igraph is provided
    2. Run Leiden clustering with component detection
    3. Calculate weighted clustering coefficients (Barrat's method)
    4. Apply broken stick model to find CC threshold per community
    5. Remove references above threshold (hubs), keep below (specific)
    6. Compact alignments and rebuild read indices
    7. Optionally export results to TSV and GraphML
    """
    
    if verbose:
        bf_nogil_logf_notime(b"CLUSTER", "Cluster-aware filtering with CC thresholding\n")
        bf_nogil_logf_notime(b"CLUSTER", "  Graph min edge weight: %u\n", graph_min_edge_weight)
        bf_nogil_logf_notime(b"CLUSTER", "  Min read count: %d\n", min_read_count)
        bf_nogil_logf_notime(b"CLUSTER", "  Leiden resolution: %.2f\n", leiden_resolution)
        bf_nogil_logf_notime(b"CLUSTER", "  Leiden max iterations: %d\n\n", leiden_max_iterations)

    if not existing_graph or not existing_graph.igraph_handle:
        bf_nogil_logf_notime(b"CLUSTER", "ERROR: No pre-built igraph provided to cluster-aware filtering\n")
        return -1

    cdef LeidenResults* leiden_results = leiden_clustering(
        existing_graph,
        ref_stats,
        leiden_resolution,
        leiden_max_iterations,
        verbose,
        thread_count,
        outlier_method,
        iforest_n_trees,
        iforest_subsample_size,
        iforest_contamination,
        iforest_random_seed,
        lof_k,
        lof_contamination,
        zscore_threshold,
        existing_graph.tsv_exact_connection_counts,
        existing_graph.tsv_co_mapping_averages,
        existing_graph.tsv_max_co_mappings,
        existing_graph.tsv_neighbor_multimap_avg,
        existing_graph.tsv_array_size
    )

    if not leiden_results:
        bf_nogil_logf_notime(b"CLUSTER", "ERROR: Leiden clustering failed\n")
        return -1

    cdef char* keep_flag = leiden_results.keep_flag

    cdef uint32_t ref_idx_loop, leiden_nnodes
    if pattern_data:
        leiden_nnodes = 0
        if leiden_results:
            leiden_nnodes = leiden_results.num_nodes

        for ref_idx_loop in range(array_size):
            if ref_idx_loop < pool.reference_count:
                if ref_idx_loop < leiden_nnodes:
                    pattern_data[ref_idx_loop].component_id = leiden_results.component_membership[ref_idx_loop]
                    pattern_data[ref_idx_loop].node_degree = leiden_results.node_degree[ref_idx_loop]
                    pattern_data[ref_idx_loop].leiden_community_id = leiden_results.community_membership[ref_idx_loop]
                    pattern_data[ref_idx_loop].leiden_community_cc = leiden_results.community_cc_values[ref_idx_loop]
                    pattern_data[ref_idx_loop].leiden_individual_cc = leiden_results.individual_cc_values[ref_idx_loop]
                    pattern_data[ref_idx_loop].leiden_cc_threshold = leiden_results.cc_threshold_values[ref_idx_loop]
                    pattern_data[ref_idx_loop].leiden_keep_flag = leiden_results.keep_flag[ref_idx_loop]
                    pattern_data[ref_idx_loop].leiden_anomaly_score = leiden_results.anomaly_scores[ref_idx_loop]
                    pattern_data[ref_idx_loop].betweenness_centrality = leiden_results.betweenness_centrality[ref_idx_loop]
                else:
                    pattern_data[ref_idx_loop].component_id = UINT32_MAX
                    pattern_data[ref_idx_loop].node_degree = 0
                    pattern_data[ref_idx_loop].leiden_community_id = UINT32_MAX
                    pattern_data[ref_idx_loop].leiden_community_cc = 0.0
                    pattern_data[ref_idx_loop].leiden_individual_cc = 0.0
                    pattern_data[ref_idx_loop].leiden_cc_threshold = 0.0
                    pattern_data[ref_idx_loop].leiden_keep_flag = 1
                    pattern_data[ref_idx_loop].leiden_anomaly_score = 0.0
                    pattern_data[ref_idx_loop].betweenness_centrality = 0.0
    cdef uint32_t ref_idx
    cdef int64_t refs_kept = 0, refs_removed = 0
    
    for ref_idx in range(array_size):
        if ref_stats[ref_idx].total_reads < min_read_count:
            continue
        
        if keep_flag[ref_idx]:
            refs_kept += 1
        else:
            refs_removed += 1
    
    if verbose:
        bf_nogil_logf_notime(b"CLUSTER", "\n=== Cluster-aware filtering results ===\n")
        bf_nogil_logf_notime(b"CLUSTER", "  References kept: %ld\n", refs_kept)
        bf_nogil_logf_notime(b"CLUSTER", "  References removed: %ld\n", refs_removed)
        bf_nogil_logf_notime(
            b"CLUSTER",
            "  Removal rate: %.1f%%\n",
            100.0 * <double>refs_removed / <double>(refs_kept + refs_removed),
        )

    cdef int64_t alignment_idx, new_alignment_idx = 0
    cdef uint32_t ref_id
    cdef int64_t surviving_alignments = 0

    for alignment_idx in range(pool.alignment_count):
        ref_id = pool.alignments[alignment_idx].reference_index
        if ref_id < array_size and keep_flag[ref_id]:
            surviving_alignments += 1

    if surviving_alignments == 0:
        bf_nogil_logf_notime(b"CLUSTER", "ERROR: No alignments survive cluster-aware filtering\n")
        free(keep_flag)
        return -1

    cdef float* new_zp_values = <float*>malloc(surviving_alignments * sizeof(float))
    if not new_zp_values:
        free(keep_flag)
        return -1
    new_alignment_idx = 0
    for alignment_idx in range(pool.alignment_count):
        ref_id = pool.alignments[alignment_idx].reference_index
        
        if ref_id < array_size and keep_flag[ref_id]:
            if new_alignment_idx != alignment_idx:
                pool.alignments[new_alignment_idx] = pool.alignments[alignment_idx]
            
            if pool.precomputed_zp_values:
                new_zp_values[new_alignment_idx] = pool.precomputed_zp_values[alignment_idx]
            
            new_alignment_idx += 1

    if pool.precomputed_zp_values:
        free(pool.precomputed_zp_values)
    pool.precomputed_zp_values = new_zp_values

    cdef int64_t alignments_removed = pool.alignment_count - new_alignment_idx
    pool.alignment_count = new_alignment_idx

    if verbose:
        bf_nogil_logf_notime(b"CLUSTER", "\n=== Alignment compaction ===\n")
        bf_nogil_logf_notime(b"CLUSTER", "  Alignments removed: %ld\n", alignments_removed)
        bf_nogil_logf_notime(b"CLUSTER", "  Alignments kept: %ld\n", pool.alignment_count)

    memset(pool.read_alignment_counts, 0, pool.unique_read_count * sizeof(uint32_t))

    cdef uint32_t current_rid
    for alignment_idx in range(pool.alignment_count):
        current_rid = pool.alignments[alignment_idx].read_index
        if current_rid < pool.unique_read_count:
            pool.read_alignment_counts[current_rid] += 1

    cdef uint64_t start_pos = 0
    cdef uint32_t rid
    for rid in range(pool.unique_read_count):
        pool.read_alignment_starts[rid] = start_pos
        start_pos += pool.read_alignment_counts[rid]
    cdef int64_t reads_surviving = 0
    for rid in range(pool.unique_read_count):
        if pool.read_alignment_counts[rid] > 0:
            reads_surviving += 1
    
    pool.final_unique_reads = reads_surviving
    
    if verbose:
        bf_nogil_logf_notime(b"CLUSTER", "\n=== Read statistics ===\n")
        bf_nogil_logf_notime(b"CLUSTER", "  Original unique reads: %u\n", pool.unique_read_count)
        bf_nogil_logf_notime(
            b"CLUSTER",
            "  Reads with surviving alignments: %ld (%.1f%%)\n",
            reads_surviving,
            100.0 * <double>reads_surviving / <double>pool.unique_read_count,
        )
        bf_nogil_logf_notime(b"CLUSTER", "=================================================================\n\n")

    if tsv_export_path and existing_graph:
        if verbose:
            bf_nogil_logf_notime(b"CLUSTER", "Writing TSV with Leiden results to: %s\n", tsv_export_path)

    if write_graph_tsv(pool, bam_header, mapping,
                          pattern_data, ref_stats,
                          existing_graph.tsv_total_reads,
                          existing_graph.tsv_multimap_reads,
                          existing_graph.tsv_alignments_per_ref,
                          existing_graph.tsv_exact_connection_counts,
                          existing_graph.tsv_co_mapping_averages,
                          existing_graph.tsv_max_co_mappings,
                          existing_graph.tsv_co_mapping_counts,
                          existing_graph.tsv_neighbor_multimap_avg,
                          existing_graph.tsv_neighbor_connections_avg,
                          existing_graph.tsv_neighbor_counts,
                          existing_graph.tsv_array_size,
                          existing_graph.tsv_dataset_median_connections,
                          existing_graph.tsv_min_read_count,
              True,
              outlier_method,
              tsv_export_path) != 0:
            bf_nogil_logf_notime(b"CLUSTER", "ERROR: Failed to write graph analysis TSV with Leiden results\n")

    if graph_export_path and existing_graph:
        if verbose:
            bf_nogil_logf_notime(b"CLUSTER", "Exporting graph to GraphML format: %s\n", graph_export_path)

        if export_graph_graphml(existing_graph, pool, bam_header, mapping,
                               pattern_data, ref_stats,
                               graph_export_path, verbose, outlier_method) != 0:
            bf_nogil_logf_notime(b"CLUSTER", "ERROR: Failed to export graph to GraphML\n")

    if leiden_results:
        if leiden_results.keep_flag:
            free(leiden_results.keep_flag)
        if leiden_results.community_membership:
            free(leiden_results.community_membership)
        if leiden_results.component_membership:
            free(leiden_results.component_membership)
        if leiden_results.node_degree:
            free(leiden_results.node_degree)
        if leiden_results.community_cc_values:
            free(leiden_results.community_cc_values)
        if leiden_results.individual_cc_values:
            free(leiden_results.individual_cc_values)
        if leiden_results.cc_threshold_values:
            free(leiden_results.cc_threshold_values)
        if leiden_results.anomaly_scores:
            free(leiden_results.anomaly_scores)
        if leiden_results.betweenness_centrality:
            free(leiden_results.betweenness_centrality)
        free(leiden_results)

    return 0
