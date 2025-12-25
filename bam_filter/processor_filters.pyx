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
from libc.stdint cimport int64_t, int32_t, uint32_t, uint64_t, uint8_t
from libc.math cimport log2, fmin, fmax, exp, sqrt, pow, INFINITY

cdef uint32_t UINT32_MAX = 0xFFFFFFFF
from bam_filter.processor cimport (
    MemoryPool,
    AlignmentCore, HierarchicalData, DamageCounts, BAMWriterAux,
)
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
    prune_low_weight_edges,
    calculate_graph_statistics,
    extract_neighbors_from_igraph
)
from bam_filter.processor_community_igraph cimport (
    community_clustering,
    CommunityResults,
)
from bam_filter.processor_igraph cimport *
from bam_filter.processor_graph cimport write_graph_tsv
from bam_filter.processor_taxonomy_filters cimport (
    TaxonomyFilterConfig,
    TaxonomyFilterStats,
    weight_anomaly_scores_by_taxonomy,
    apply_taxonomy_informed_filtering
)
from bam_filter.processor_graph_export cimport export_graph_graphml
from bam_filter.processor_mapping cimport ReferenceMapping
from bam_filter.taxonomy_db cimport TaxonomyDB
from bam_filter.processor_tiered_filters cimport apply_tiered_filtering
from bam_filter.processor_graph_taxonomy cimport (
    detect_taxonomy_anomalies,
    TaxonomyGraphConfig
)
from bam_filter.processor_network_qc cimport (
    NetworkQCConfig,
    NetworkQCReference,
    compute_reference_qc_metrics,
    compute_tax_ambiguity_flags,
)

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
    Supports both legacy alignments array and split array storage. With split arrays,
    ANI filtering was already done during population so all alignments pass.
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
    cdef bint use_split = (pool.alignment_cores != NULL)

    precomp = create_precomputed_weights(pool.reference_count)
    if not precomp:
        return -1
    update_precomputed_weights(precomp, reference_weights)

    bf_nogil_logf_notime(
        b"FILTER",
        "probability_filter: start pmd_output=%s total_alignments=%lld split=%d",
        b"enabled" if pool.pmd_enabled_for_output else b"disabled",
        <long long>pool.alignment_count,
        <int>use_split,
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

    # Phase 1: Count survivors per read
    for read_idx in prange(pool.unique_read_count, nogil=True, schedule='static', num_threads=config.thread_count):
        alignment_count = pool.read_alignment_counts[read_idx]
        if alignment_count == 0:
            continue
        start_pos = pool.read_alignment_starts[read_idx]
        end_pos = start_pos + alignment_count
        if alignment_count == 1:
            read_max_probs[read_idx] = 1.0
            # With split arrays, all alignments already passed ANI filter
            if 1.0 >= min_threshold and (fraction_threshold == 0.0 or 1.0 >= fraction_threshold * 1.0):
                survivors_per_read[read_idx] = 1
            continue

        log_norm = NEG_INF
        for ai in range(start_pos, end_pos):
            if use_split:
                ref_idx = pool.alignment_cores[ai].reference_index
                alignment_score = pool.alignment_cores[ai].alignment_score
            else:
                ref_idx = pool.alignments[ai].reference_index
                alignment_score = pool.alignments[ai].alignment_score
            if ref_idx < pool.reference_count:
                log_lik = <double>alignment_score
                log_weighted = precomp.log_weights[ref_idx] + log_lik
                log_norm = stable_log_sum_exp(log_norm, log_weighted)

        if log_norm == NEG_INF:
            uniform_zp = 1.0 / alignment_count if alignment_count > 0 else 0.0
            read_max_probs[read_idx] = uniform_zp
            if uniform_zp >= min_threshold and (fraction_threshold == 0.0 or uniform_zp >= fraction_threshold * uniform_zp):
                # With split arrays, all alignments already passed ANI filter
                survivors_per_read[read_idx] = alignment_count
            continue

        p = 0.0
        for ai in range(start_pos, end_pos):
            if use_split:
                ref_idx = pool.alignment_cores[ai].reference_index
                alignment_score = pool.alignment_cores[ai].alignment_score
            else:
                ref_idx = pool.alignments[ai].reference_index
                alignment_score = pool.alignments[ai].alignment_score
            if ref_idx < pool.reference_count:
                log_lik = <double>alignment_score
                log_weighted = precomp.log_weights[ref_idx] + log_lik
                posterior = exp(log_weighted - log_norm)
                if posterior > p:
                    p = posterior
        read_max_probs[read_idx] = <float>p
        thr = fraction_threshold * p
        for ai in range(start_pos, end_pos):
            if use_split:
                ref_idx = pool.alignment_cores[ai].reference_index
                alignment_score = pool.alignment_cores[ai].alignment_score
            else:
                ref_idx = pool.alignments[ai].reference_index
                alignment_score = pool.alignments[ai].alignment_score
            if ref_idx < pool.reference_count:
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

    # Phase 2: Compact alignments
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
                    if use_split:
                        pool.alignment_cores[write_idx] = pool.alignment_cores[start_pos]
                        pool.read_indices[write_idx] = pool.read_indices[start_pos]
                        if pool.hierarchical != NULL:
                            pool.hierarchical[write_idx] = pool.hierarchical[start_pos]
                        if pool.damage_counts != NULL:
                            pool.damage_counts[write_idx] = pool.damage_counts[start_pos]
                        if pool.bam_aux != NULL:
                            pool.bam_aux[write_idx] = pool.bam_aux[start_pos]
                    else:
                        pool.alignments[write_idx] = pool.alignments[start_pos]
                pool.precomputed_zp_values[write_idx] = 1.0
                write_idx += 1
            continue

        log_norm = NEG_INF
        for ai in range(start_pos, end_pos):
            if use_split:
                ref_idx = pool.alignment_cores[ai].reference_index
                alignment_score = pool.alignment_cores[ai].alignment_score
            else:
                ref_idx = pool.alignments[ai].reference_index
                alignment_score = pool.alignments[ai].alignment_score
            if ref_idx < pool.reference_count:
                log_weighted = precomp.log_weights[ref_idx] + alignment_score
                log_norm = stable_log_sum_exp(log_norm, log_weighted)

        if log_norm == NEG_INF:
            uniform_zp = 1.0 / alignment_count if alignment_count > 0 else 0.0
            if uniform_zp >= min_threshold and (fraction_threshold == 0.0 or uniform_zp >= fraction_threshold * uniform_zp):
                for ai in range(start_pos, end_pos):
                    if write_idx != ai:
                        if use_split:
                            pool.alignment_cores[write_idx] = pool.alignment_cores[ai]
                            pool.read_indices[write_idx] = pool.read_indices[ai]
                            if pool.hierarchical != NULL:
                                pool.hierarchical[write_idx] = pool.hierarchical[ai]
                            if pool.damage_counts != NULL:
                                pool.damage_counts[write_idx] = pool.damage_counts[ai]
                            if pool.bam_aux != NULL:
                                pool.bam_aux[write_idx] = pool.bam_aux[ai]
                        else:
                            pool.alignments[write_idx] = pool.alignments[ai]
                    pool.precomputed_zp_values[write_idx] = uniform_zp
                    write_idx += 1
            continue

        thr = (<double>fraction_threshold) * (<double>read_max_probs[read_idx])
        for ai in range(start_pos, end_pos):
            if use_split:
                ref_idx = pool.alignment_cores[ai].reference_index
                alignment_score = pool.alignment_cores[ai].alignment_score
            else:
                ref_idx = pool.alignments[ai].reference_index
                alignment_score = pool.alignments[ai].alignment_score
            if ref_idx < pool.reference_count:
                log_weighted = precomp.log_weights[ref_idx] + alignment_score
                posterior = exp(log_weighted - log_norm)
                keep_alignment = (posterior >= min_threshold) and (fraction_threshold == 0.0 or posterior >= thr)
                if keep_alignment:
                    if write_idx != ai:
                        if use_split:
                            pool.alignment_cores[write_idx] = pool.alignment_cores[ai]
                            pool.read_indices[write_idx] = pool.read_indices[ai]
                            if pool.hierarchical != NULL:
                                pool.hierarchical[write_idx] = pool.hierarchical[ai]
                            if pool.damage_counts != NULL:
                                pool.damage_counts[write_idx] = pool.damage_counts[ai]
                            if pool.bam_aux != NULL:
                                pool.bam_aux[write_idx] = pool.bam_aux[ai]
                        else:
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

    # Rebuild read index structures
    memset(pool.read_alignment_counts, 0, pool.unique_read_count * sizeof(uint32_t))
    for ai in range(pool.alignment_count):
        if use_split:
            rid = pool.read_indices[ai]
        else:
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
    Supports both legacy alignments array and split array storage.
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
    cdef bint use_split = (pool.alignment_cores != NULL)

    if not reference_counts or not reference_keep_flag:
        if reference_counts: free(reference_counts)
        if reference_keep_flag: free(reference_keep_flag)
        return -1

    for alignment_idx in range(pool.alignment_count):
        if use_split:
            ref_id = pool.alignment_cores[alignment_idx].reference_index
        else:
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
        if use_split:
            ref_id = pool.alignment_cores[alignment_idx].reference_index
        else:
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
        if use_split:
            ref_id = pool.alignment_cores[alignment_idx].reference_index
        else:
            ref_id = pool.alignments[alignment_idx].reference_index

        if (ref_id < pool.reference_count and reference_keep_flag[ref_id]):
            if new_alignment_idx != alignment_idx:
                if use_split:
                    pool.alignment_cores[new_alignment_idx] = pool.alignment_cores[alignment_idx]
                    pool.read_indices[new_alignment_idx] = pool.read_indices[alignment_idx]
                    if pool.hierarchical != NULL:
                        pool.hierarchical[new_alignment_idx] = pool.hierarchical[alignment_idx]
                    if pool.damage_counts != NULL:
                        pool.damage_counts[new_alignment_idx] = pool.damage_counts[alignment_idx]
                    if pool.bam_aux != NULL:
                        pool.bam_aux[new_alignment_idx] = pool.bam_aux[alignment_idx]
                else:
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

    if alignments_removed > 0:
        bf_nogil_logf_notime(b"FILTER-PMD", "Rebuilding read indices...\n")

        memset(pool.read_alignment_counts, 0, pool.unique_read_count * sizeof(uint32_t))
        for alignment_idx in range(pool.alignment_count):
            if use_split:
                current_rid = pool.read_indices[alignment_idx]
            else:
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
                                       bint verbose,
                                       double community_resolution,
                                       int community_max_iterations,
                                       uint32_t graph_min_edge_weight,
                                       int32_t thread_count,
                                       int outlier_method,
                                       WeightedGraph* existing_graph,
                                       sam_hdr_t* bam_header,
                                       ReferenceMapping* mapping,
                                       const char* tsv_export_path,
                                       const char* graph_export_path,
                                       TaxonomyFilterConfig* taxonomy_filter_config,
                                       TaxonomyFilterStats* taxonomy_stats_out,
                                       TaxonomyDB* taxonomy_db,
                                       float betweenness_threshold,
                                       float cc_threshold,
                                       uint32_t hub_degree_threshold,
                                       bint strict_mode,
                                       bint remove_cross_domain_edges,
                                       bint flag_misannotations,
                                       NetworkQCConfig* network_qc_config,
                                       uint8_t tax_ambiguity_removal_level) except -1 nogil:
    """Apply three-tier filtering that combines clustering, topology metrics, and taxonomy checks.

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
    verbose : bint
        Enable verbose logging
    community_resolution : double
        Community resolution parameter
    community_max_iterations : int
        Maximum Community iterations
    graph_min_edge_weight : uint32_t
        Minimum edge weight for graph construction
    thread_count : int32_t
        Thread count for parallel processing
    outlier_method : int
        Outlier detection method (MAD recommended)
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
    taxonomy_filter_config : TaxonomyFilterConfig*
        Taxonomy filtering configuration
    taxonomy_stats_out : TaxonomyFilterStats*
        Taxonomy filtering statistics (output)
    taxonomy_db : TaxonomyDB*
        Taxonomy database for coherence checking
    betweenness_threshold : float
        Threshold for bridge detection (e.g., 0.01)
    cc_threshold : float
        Threshold for hub detection (e.g., 0.3)
    hub_degree_threshold : uint32_t
        Minimum degree for hub classification (e.g., 5)
    strict_mode : bint
        Enable strict filtering mode

    Returns
    -------
    int
        0 on success, -1 on error

    Notes
    -----
    Workflow:
    1. Run community clustering with component detection.
    2. Compute betweenness and clustering coefficients.
    3. Apply taxonomy-informed filtering when enabled.
    4. Execute the tiered decision pipeline (structural role, coherence, integrated decision).
    5. Compact alignments, rebuild indices, and emit optional exports.
    """
    
    if verbose:
        bf_nogil_logf_notime(b"CLUSTER", "Cluster-aware filtering with CC thresholding\n")
        bf_nogil_logf_notime(b"CLUSTER", "  Graph min edge weight: %u\n", graph_min_edge_weight)
        bf_nogil_logf_notime(b"CLUSTER", "  Min read count: %d\n", min_read_count)
        # Show algorithm-specific parameters
        # Note: Resolution and max_iterations only apply to Community, not to LPA
        # LPA will be used for graphs >= 10k nodes
        if existing_graph.num_nodes >= 10000:
            bf_nogil_logf_notime(b"CLUSTER", "  Algorithm: Label Propagation (LPA) - fast O(m) for large graphs\n")
        else:
            bf_nogil_logf_notime(b"CLUSTER", "  Algorithm: Leiden (community detection)\n")
            bf_nogil_logf_notime(b"CLUSTER", "  Community resolution: %.2f\n", community_resolution)
            bf_nogil_logf_notime(b"CLUSTER", "  Community max iterations: %d\n", community_max_iterations)
        bf_nogil_logf_notime(b"CLUSTER", "\n")

    if not existing_graph or not existing_graph.igraph_handle:
        bf_nogil_logf_notime(b"CLUSTER", "ERROR: No pre-built igraph provided to cluster-aware filtering\n")
        return -1

    # Level 2 taxonomy filtering: Weight anomaly scores BEFORE Community clustering
    # This array will be allocated and filled by community_clustering, but we need to
    # create a placeholder anomaly_scores array that will be passed to weight_anomaly_scores_by_taxonomy
    # Note: The actual anomaly score computation happens inside community_clustering, so we'll apply
    # taxonomy weighting AFTER community_clustering returns the anomaly scores. This is handled by passing
    # the taxonomy_filter_config to community_clustering.

    if verbose:
        bf_nogil_logf_notime(b"CLUSTER", "Calling community_clustering with %u nodes...\n", existing_graph.num_nodes)

    cdef CommunityResults* community_results = community_clustering(
        existing_graph,
        ref_stats,
        community_resolution,
        community_max_iterations,
        verbose,
        thread_count,
        outlier_method,
        existing_graph.tsv_exact_connection_counts,
        existing_graph.tsv_co_mapping_averages,
        existing_graph.tsv_max_co_mappings,
        existing_graph.tsv_neighbor_multimap_avg,
        existing_graph.tsv_array_size
    )

    if not community_results:
        bf_nogil_logf_notime(b"CLUSTER", "ERROR: Community clustering failed\n")
        return -1


    cdef char* keep_flag = community_results.keep_flag

    # Detect taxonomy anomalies using neighbor connectivity (Phase 6b)
    # This requires that taxonomy IDs were already populated in Phase 5b
    cdef TaxonomyGraphConfig tax_config
    cdef uint32_t** neighbor_lists = NULL
    cdef uint32_t* neighbor_counts = NULL
    cdef uint32_t tax_ref_idx
    cdef int extract_result = 0
    cdef bint neighbors_from_igraph = False


    if taxonomy_filter_config != NULL and taxonomy_filter_config.enabled and taxonomy_db != NULL:
        if verbose:
            bf_nogil_logf_notime(b"CLUSTER", "Detecting taxonomy anomalies using graph connectivity...\n")

        # Configure taxonomy anomaly detection
        # Use reasonable defaults for detection thresholds (not filtering thresholds)
        tax_config.enabled = True
        tax_config.min_rank_id_for_comparison = 6  # genus level
        tax_config.cross_domain_threshold = 0.10  # Flag if >10% of neighbors are cross-domain
        tax_config.kingdom_mismatch_threshold = 0.25  # Flag if >25% are kingdom mismatches
        tax_config.genus_mismatch_threshold = 0.50  # Flag if >50% are genus+ mismatches

        # Extract neighbor lists from existing_graph
        # Two cases: nodes array exists, or extract from igraph
        if existing_graph:
            if existing_graph.nodes:
                # Case 1: Full WeightedGraph with nodes array
                neighbor_lists = <uint32_t**>malloc(array_size * sizeof(uint32_t*))
                neighbor_counts = <uint32_t*>malloc(array_size * sizeof(uint32_t))
                if neighbor_lists and neighbor_counts:
                    for tax_ref_idx in range(array_size):
                        if tax_ref_idx >= existing_graph.num_nodes:
                            neighbor_lists[tax_ref_idx] = NULL
                            neighbor_counts[tax_ref_idx] = 0
                        else:
                            neighbor_lists[tax_ref_idx] = existing_graph.nodes[tax_ref_idx].neighbors
                            neighbor_counts[tax_ref_idx] = existing_graph.nodes[tax_ref_idx].degree
                else:
                    bf_nogil_logf_notime(b"CLUSTER", "WARNING: Failed to allocate neighbor arrays\n")
            elif existing_graph.igraph_handle:
                # Case 2: Extract from igraph (clustering mode)
                extract_result = extract_neighbors_from_igraph(
                    existing_graph, array_size,
                    &neighbor_lists, &neighbor_counts, verbose)
                if extract_result == 0:
                    neighbors_from_igraph = True
                else:
                    bf_nogil_logf_notime(b"CLUSTER", "WARNING: Failed to extract neighbors from igraph\n")

            # Run taxonomy anomaly detection if we have neighbor data
            if neighbor_lists and neighbor_counts:
                detect_taxonomy_anomalies(
                    pattern_data,
                    array_size,
                    pool,
                    taxonomy_db,
                    &tax_config,
                    neighbor_lists,
                    neighbor_counts,
                    <void*>existing_graph,
                    verbose
                )

                # NOTE: neighbor_lists and neighbor_counts cleanup moved to AFTER
                # apply_tiered_filtering so they can be used for edge removal

                if verbose:
                    bf_nogil_logf_notime(b"CLUSTER", "Taxonomy anomaly detection complete\n")
        else:
            bf_nogil_logf_notime(b"CLUSTER", "WARNING: No graph available for taxonomy anomaly detection\n")

    # Apply taxonomy-informed filtering (Levels 1 & 3) AFTER Community clustering
    # Note: Level 2 (weighted anomaly scores) is skipped in this integration because
    # anomaly scores are computed inside community_clustering. To fully integrate Level 2,
    # we would need to modify community_clustering to accept taxonomy_filter_config.
    cdef TaxonomyFilterStats taxonomy_stats
    taxonomy_stats.strict_removed = 0
    taxonomy_stats.weighted_count = 0
    taxonomy_stats.second_chance_restored = 0

    if taxonomy_filter_config != NULL and taxonomy_filter_config.enabled:
        if verbose:
            bf_nogil_logf_notime(b"CLUSTER", "Applying taxonomy-informed filtering...\n")

        # Apply strict filtering (Level 1) and second-chance validation (Level 3)
        taxonomy_stats = apply_taxonomy_informed_filtering(
            pattern_data,
            array_size,
            keep_flag,
            community_results.anomaly_scores,
            existing_graph.tsv_exact_connection_counts,
            taxonomy_filter_config,
            verbose
        )

    cdef uint32_t ref_idx_loop, community_nnodes
    if pattern_data:
        community_nnodes = 0
        if community_results:
            community_nnodes = community_results.num_nodes

        for ref_idx_loop in range(array_size):
            if ref_idx_loop < pool.reference_count:
                if ref_idx_loop < community_nnodes:
                    pattern_data[ref_idx_loop].component_id = community_results.component_membership[ref_idx_loop]
                    pattern_data[ref_idx_loop].node_degree = community_results.node_degree[ref_idx_loop]
                    # Store original_degree from graph (before pruning) for pruned isolated node detection
                    if existing_graph and existing_graph.nodes and ref_idx_loop < existing_graph.num_nodes:
                        pattern_data[ref_idx_loop].original_degree = existing_graph.nodes[ref_idx_loop].original_degree
                    else:
                        pattern_data[ref_idx_loop].original_degree = pattern_data[ref_idx_loop].node_degree
                    pattern_data[ref_idx_loop].community_id = community_results.community_membership[ref_idx_loop]
                    pattern_data[ref_idx_loop].community_cc = community_results.community_cc_values[ref_idx_loop]
                    pattern_data[ref_idx_loop].community_individual_cc = community_results.individual_cc_values[ref_idx_loop]
                    pattern_data[ref_idx_loop].community_cc_threshold = community_results.cc_threshold_values[ref_idx_loop]
                    pattern_data[ref_idx_loop].community_keep_flag = community_results.keep_flag[ref_idx_loop]
                    pattern_data[ref_idx_loop].community_anomaly_score = community_results.anomaly_scores[ref_idx_loop]
                    pattern_data[ref_idx_loop].betweenness_centrality = community_results.betweenness_centrality[ref_idx_loop]
                    pattern_data[ref_idx_loop].num_neighbor_communities = community_results.num_neighbor_communities[ref_idx_loop]
                else:
                    pattern_data[ref_idx_loop].component_id = UINT32_MAX
                    pattern_data[ref_idx_loop].node_degree = 0
                    pattern_data[ref_idx_loop].original_degree = 0
                    pattern_data[ref_idx_loop].community_id = UINT32_MAX
                    pattern_data[ref_idx_loop].community_cc = 0.0
                    pattern_data[ref_idx_loop].community_individual_cc = 0.0
                    pattern_data[ref_idx_loop].community_cc_threshold = 0.0
                    pattern_data[ref_idx_loop].community_keep_flag = 1
                    pattern_data[ref_idx_loop].community_anomaly_score = 0.0
                    pattern_data[ref_idx_loop].betweenness_centrality = 0.0
                    pattern_data[ref_idx_loop].num_neighbor_communities = 0


    # Apply three-tier filtering (NEW - includes bridge detection!)
    if verbose:
        bf_nogil_logf_notime(b"CLUSTER", "\\nApplying three-tier filtering...\\n")

    # Use thresholds passed as function parameters
    # Pass neighbor data for edge removal (re-extract from graph if needed)
    cdef char* alignment_keep_flags = NULL
    if apply_tiered_filtering(
        pattern_data,
        ref_stats,
        array_size,
        keep_flag,
        taxonomy_db,
        betweenness_threshold,
        cc_threshold,
        hub_degree_threshold,
        strict_mode,
        verbose,
        <void*>pool,
        neighbor_lists,
        neighbor_counts,
        remove_cross_domain_edges,  # enable_edge_removal from args
        flag_misannotations,
        &alignment_keep_flags  # Get alignment flags for combined compaction
    ) != 0:
        bf_nogil_logf_notime(b"CLUSTER", "ERROR: Tiered filtering failed\\n")
        # Continue anyway, don't fail completely

    # =========================================================================
    # Network QC Filtering (Taxonomic Ambiguity Detection)
    # =========================================================================
    # Computes neighbor taxonomy entropy per node and flags references with
    # taxonomically diverse neighbors. HUB references are skipped since they
    # are expected to have high entropy by nature (many connections).
    cdef NetworkQCReference* network_qc_metrics = NULL
    cdef uint32_t tax_ambiguity_removed = 0
    cdef uint32_t qc_idx
    cdef int32_t* community_ids_for_qc = NULL
    cdef int32_t* taxonomy_ids_for_qc = NULL
    cdef char* structural_roles_for_qc = NULL
    cdef uint32_t max_community_id = 0
    cdef uint32_t nodes_with_neighbors = 0
    cdef uint32_t nodes_processed = 0
    cdef uint32_t flag_counts[4]
    cdef uint32_t fl

    if network_qc_config != NULL and tax_ambiguity_removal_level > 0:
        if verbose:
            bf_nogil_logf_notime(b"NETQC", "Applying network QC filtering (removal_level=%d)...\n", tax_ambiguity_removal_level)

        # Allocate network QC metrics array
        network_qc_metrics = <NetworkQCReference*>calloc(array_size, sizeof(NetworkQCReference))
        if network_qc_metrics == NULL:
            bf_nogil_logf_notime(b"NETQC", "ERROR: Failed to allocate network QC metrics\\n")
        else:
            # Build community IDs, taxonomy IDs, and structural roles arrays
            community_ids_for_qc = <int32_t*>calloc(array_size, sizeof(int32_t))
            taxonomy_ids_for_qc = <int32_t*>calloc(array_size, sizeof(int32_t))
            structural_roles_for_qc = <char*>calloc(array_size, sizeof(char))

            if community_ids_for_qc != NULL and taxonomy_ids_for_qc != NULL and structural_roles_for_qc != NULL:
                # Populate community and taxonomy IDs
                # Note: UINT32_MAX is used as sentinel for isolated nodes
                for qc_idx in range(array_size):
                    if community_results != NULL and qc_idx < community_results.num_nodes:
                        if community_results.community_membership[qc_idx] == UINT32_MAX:
                            # Isolated node - use -1 for int32 representation
                            community_ids_for_qc[qc_idx] = -1
                        else:
                            community_ids_for_qc[qc_idx] = <int32_t>community_results.community_membership[qc_idx]
                            # Track max valid community ID (skip sentinel values)
                            if community_results.community_membership[qc_idx] > max_community_id:
                                max_community_id = community_results.community_membership[qc_idx]
                    else:
                        community_ids_for_qc[qc_idx] = -1  # No community

                    # Get taxonomy ID at domain level for entropy computation
                    # Using domain-level taxids prevents over-flagging due to species-level diversity
                    if pattern_data != NULL:
                        # Use domain_taxid for entropy (species-level taxids cause over-filtering)
                        if pattern_data[qc_idx].domain_taxid > 0:
                            taxonomy_ids_for_qc[qc_idx] = pattern_data[qc_idx].domain_taxid
                        else:
                            # Fallback to species-level taxid if domain not found
                            taxonomy_ids_for_qc[qc_idx] = pattern_data[qc_idx].taxid
                        structural_roles_for_qc[qc_idx] = pattern_data[qc_idx].structural_role
                    else:
                        taxonomy_ids_for_qc[qc_idx] = -1
                        structural_roles_for_qc[qc_idx] = 0  # PERIPHERAL

                # Compute per-reference QC metrics
                # Note: This is a simplified version that uses existing neighbor data
                nodes_with_neighbors = 0
                nodes_processed = 0
                for qc_idx in range(array_size):
                    if neighbor_lists != NULL and neighbor_counts != NULL and neighbor_lists[qc_idx] != NULL:
                        nodes_with_neighbors += 1
                        if neighbor_counts[qc_idx] > 0:
                            nodes_processed += 1
                            compute_reference_qc_metrics(
                                qc_idx,
                                &network_qc_metrics[qc_idx],
                                community_ids_for_qc,
                                taxonomy_ids_for_qc,
                                neighbor_lists[qc_idx],
                                NULL,  # neighbor_weights (use degree)
                                neighbor_counts[qc_idx],
                                NULL,  # community_strength_sums
                                NULL,  # community_strength_sq
                                NULL,  # community_counts
                                max_community_id + 1,  # n_communities
                                community_ids_for_qc[qc_idx],  # ref_community
                                network_qc_config,
                            )

                # Compute taxonomic ambiguity flags (skips HUB references)
                compute_tax_ambiguity_flags(network_qc_metrics, array_size, structural_roles_for_qc, network_qc_config)

                # Count tax ambiguity flags at each level for logging
                for fl in range(4):
                    flag_counts[fl] = 0
                for qc_idx in range(array_size):
                    if network_qc_metrics[qc_idx].tax_ambiguity_flag < 4:
                        flag_counts[network_qc_metrics[qc_idx].tax_ambiguity_flag] += 1

                if verbose:
                    bf_nogil_logf_notime(b"NETQC", "qc_analyzed=%u communities=%u flags: clean=%u biased=%u mixed=%u highly_mixed=%u\n",
                                         nodes_processed, max_community_id + 1,
                                         flag_counts[0], flag_counts[1], flag_counts[2], flag_counts[3])

                # Apply taxonomic ambiguity filtering and copy metrics to pattern_data for TSV export
                for qc_idx in range(array_size):
                    # Copy QC metrics to pattern_data for export
                    if pattern_data != NULL:
                        pattern_data[qc_idx].neighbor_tax_entropy = network_qc_metrics[qc_idx].neighbor_tax_entropy
                        pattern_data[qc_idx].tax_ambiguity_flag = network_qc_metrics[qc_idx].tax_ambiguity_flag

                    # Apply filtering based on tax ambiguity flag
                    if network_qc_metrics[qc_idx].tax_ambiguity_flag >= tax_ambiguity_removal_level:
                        if keep_flag[qc_idx]:  # Only count if not already removed
                            keep_flag[qc_idx] = 0
                            tax_ambiguity_removed += 1

                bf_nogil_logf_notime(b"NETQC", "Network QC filtering: removed %u taxonomically ambiguous references\n", tax_ambiguity_removed)

            # Cleanup
            if community_ids_for_qc != NULL:
                free(community_ids_for_qc)
            if taxonomy_ids_for_qc != NULL:
                free(taxonomy_ids_for_qc)
            if structural_roles_for_qc != NULL:
                free(structural_roles_for_qc)
            free(network_qc_metrics)

    # Cleanup neighbor lists (after edge removal is complete)
    if neighbor_lists != NULL and neighbor_counts != NULL:
        if neighbors_from_igraph:
            # Deep free - we allocated individual arrays
            for tax_ref_idx in range(array_size):
                if neighbor_lists[tax_ref_idx]:
                    free(neighbor_lists[tax_ref_idx])
        # Free outer arrays
        free(neighbor_lists)
        free(neighbor_counts)
        neighbor_lists = NULL
        neighbor_counts = NULL

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
    cdef int64_t edge_removed_count = 0
    cdef int64_t ref_removed_count = 0
    cdef bint use_split = (pool.alignment_cores != NULL)

    # Count surviving alignments (must pass BOTH filters)
    for alignment_idx in range(pool.alignment_count):
        if use_split:
            ref_id = pool.alignment_cores[alignment_idx].reference_index
        else:
            ref_id = pool.alignments[alignment_idx].reference_index

        # Check reference filter
        if ref_id >= array_size or not keep_flag[ref_id]:
            ref_removed_count += 1
            continue

        # Check edge removal filter (if enabled)
        if alignment_keep_flags != NULL and not alignment_keep_flags[alignment_idx]:
            edge_removed_count += 1
            continue

        surviving_alignments += 1

    if surviving_alignments == 0:
        bf_nogil_logf_notime(b"CLUSTER", "ERROR: No alignments survive cluster-aware filtering\n")
        if alignment_keep_flags: free(alignment_keep_flags)
        free(keep_flag)
        return -1

    cdef float* new_zp_values = <float*>malloc(surviving_alignments * sizeof(float))
    if not new_zp_values:
        if alignment_keep_flags: free(alignment_keep_flags)
        free(keep_flag)
        return -1

    # COMBINED COMPACTION: Apply both filters in single pass
    new_alignment_idx = 0
    for alignment_idx in range(pool.alignment_count):
        if use_split:
            ref_id = pool.alignment_cores[alignment_idx].reference_index
        else:
            ref_id = pool.alignments[alignment_idx].reference_index

        # Check reference filter
        if ref_id >= array_size or not keep_flag[ref_id]:
            continue

        # Check edge removal filter (if enabled)
        if alignment_keep_flags != NULL and not alignment_keep_flags[alignment_idx]:
            continue

        # Alignment survives both filters - keep it
        if new_alignment_idx != alignment_idx:
            if use_split:
                # Copy all split arrays
                pool.alignment_cores[new_alignment_idx] = pool.alignment_cores[alignment_idx]
                pool.read_indices[new_alignment_idx] = pool.read_indices[alignment_idx]
                if pool.hierarchical != NULL:
                    pool.hierarchical[new_alignment_idx] = pool.hierarchical[alignment_idx]
                if pool.damage_counts != NULL:
                    pool.damage_counts[new_alignment_idx] = pool.damage_counts[alignment_idx]
                if pool.bam_aux != NULL:
                    pool.bam_aux[new_alignment_idx] = pool.bam_aux[alignment_idx]
            else:
                pool.alignments[new_alignment_idx] = pool.alignments[alignment_idx]

        if pool.precomputed_zp_values:
            new_zp_values[new_alignment_idx] = pool.precomputed_zp_values[alignment_idx]

        new_alignment_idx += 1

    # Free alignment_keep_flags now that we're done with it
    if alignment_keep_flags:
        free(alignment_keep_flags)

    if pool.precomputed_zp_values:
        free(pool.precomputed_zp_values)
    pool.precomputed_zp_values = new_zp_values

    cdef int64_t alignments_removed = pool.alignment_count - new_alignment_idx
    pool.alignment_count = new_alignment_idx

    if verbose:
        bf_nogil_logf_notime(b"CLUSTER", "\n=== COMBINED ALIGNMENT COMPACTION ===\n")
        bf_nogil_logf_notime(b"CLUSTER", "  Removed by reference filtering: %ld\n", ref_removed_count)
        bf_nogil_logf_notime(b"CLUSTER", "  Removed by edge filtering:      %ld\n", edge_removed_count)
        bf_nogil_logf_notime(b"CLUSTER", "  Total alignments removed:       %ld\n", alignments_removed)
        bf_nogil_logf_notime(b"CLUSTER", "  Alignments kept:                %ld\n", pool.alignment_count)
        bf_nogil_logf_notime(b"CLUSTER", "  Read loss: %.1f%%\n",
            100.0 * <double>alignments_removed / <double>(pool.alignment_count + alignments_removed)
        )

    memset(pool.read_alignment_counts, 0, pool.unique_read_count * sizeof(uint32_t))

    cdef uint32_t current_rid
    for alignment_idx in range(pool.alignment_count):
        if use_split:
            current_rid = pool.read_indices[alignment_idx]
        else:
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
            bf_nogil_logf_notime(b"CLUSTER", "Writing TSV with Community results to: %s\n", tsv_export_path)

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
              tsv_export_path,
              taxonomy_db) != 0:
            bf_nogil_logf_notime(b"CLUSTER", "ERROR: Failed to write graph analysis TSV with Community results\n")

    if graph_export_path and existing_graph:
        if verbose:
            bf_nogil_logf_notime(b"CLUSTER", "Exporting graph to GraphML format: %s\n", graph_export_path)

        if export_graph_graphml(existing_graph, pool, bam_header, mapping,
                               pattern_data, ref_stats,
                               graph_export_path, verbose, outlier_method, taxonomy_db) != 0:
            bf_nogil_logf_notime(b"CLUSTER", "ERROR: Failed to export graph to GraphML\n")

    if community_results:
        if community_results.keep_flag:
            free(community_results.keep_flag)
        if community_results.community_membership:
            free(community_results.community_membership)
        if community_results.component_membership:
            free(community_results.component_membership)
        if community_results.node_degree:
            free(community_results.node_degree)
        if community_results.community_cc_values:
            free(community_results.community_cc_values)
        if community_results.individual_cc_values:
            free(community_results.individual_cc_values)
        if community_results.cc_threshold_values:
            free(community_results.cc_threshold_values)
        if community_results.anomaly_scores:
            free(community_results.anomaly_scores)
        if community_results.betweenness_centrality:
            free(community_results.betweenness_centrality)
        free(community_results)

    # Populate output parameter with taxonomy filtering statistics
    if taxonomy_stats_out != NULL:
        taxonomy_stats_out[0] = taxonomy_stats

    return 0
