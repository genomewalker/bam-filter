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
Fast Leiden clustering using igraph C library with clustering coefficient filtering.

WORKFLOW:
1. Identify connected components in pre-built igraph
2. Skip singletons and pairs (auto-keep, no filtering needed)
3. For each multi-node component, run Leiden clustering (PARALLEL)
4. For each Leiden community, calculate weighted clustering coefficients (Barrat)
5. Apply adaptive percentile threshold to detect low-CC outliers
6. KEEP references with CC >= threshold (cohesive/informative references)
7. REMOVE references with CC < threshold (hubs/stars, likely contamination)

FILTERING LOGIC:
- Uses adaptive percentile thresholds (5th-10th percentile) based on distribution shape
- Default: 5th percentile (removes only bottom 5%, very conservative)
- High CC = reference's neighbors are connected to each other (clique/triangle pattern) → KEEP
- Low CC = reference connects unrelated references (hub/star pattern) → REMOVE as contamination
- More conservative than edge weight filtering (which uses broken stick on skewed distributions)
- Threshold represents the minimum "cohesiveness" required to be considered informative
"""

from libc.stdlib cimport malloc, calloc, free, realloc, qsort
from libc.time cimport clock, CLOCKS_PER_SEC
from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t, uint8_t
from libc.string cimport memset
from cython.parallel cimport prange, threadid

from bam_filter.processor_graph cimport ReferenceStats
from bam_filter.processor_graph_ops cimport WeightedGraph
from bam_filter.processor_igraph cimport *

# Import new multi-metric outlier detection module
from bam_filter.processor_outlier_detection cimport (
    ReferenceFeatures,
    OutlierDetectionResult,
    OutlierMethod,
    OUTLIER_MAD,
    OUTLIER_IQR,
    OUTLIER_IFOREST,
    OUTLIER_LOF,
    OUTLIER_ZSCORE,
    detect_outliers_multivariate_c,
    detect_outliers_iforest_c,
    detect_outliers_lof_c,
    free_outlier_result,
    extract_features_from_graph_metrics,
)

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil
    int bf_should_log(int level) nogil
    enum:
        BF_LOG_LEVEL_TRACE


# ==============================================================================
# HELPER: Access igraph vector elements
# ==============================================================================

cdef inline igraph_integer_t get_vector_int_element(igraph_vector_int_t *v, igraph_integer_t i) nogil:
    """Access element from igraph_vector_int_t."""
    return (<igraph_integer_t**>v)[0][i]


cdef inline igraph_real_t get_vector_element(igraph_vector_t *v, igraph_integer_t i) nogil:
    """Access element from igraph_vector_t."""
    return (<igraph_real_t**>v)[0][i]


# ==============================================================================
# WEIGHTED CLUSTERING COEFFICIENTS (Barrat's Method)
# ==============================================================================


cdef int calculate_weighted_clustering_coefficients_for_vids(
    igraph_t* graph,
    igraph_vector_t* edge_weights,
    igraph_vector_int_t* vids_vector,
    float* clustering_coefficients,
    uint32_t n_members
) except -1 nogil:
    """
    Calculate Barrat weighted clustering coefficients for a subset of vertices
    specified by an igraph_vector_int_t (vertex ids). This uses igraph_vs_vector
    to create a vertex selector and calls igraph_transitivity_barrat on the full
    graph with that selector. The result is copied into `clustering_coefficients`
    (length n_members).
    """
    cdef int ret
    cdef igraph_vector_t result
    cdef igraph_vs_t vs
    cdef uint32_t i

    ret = igraph_vector_init(&result, <igraph_integer_t>n_members)
    if ret != 0:
        return -1

    ret = igraph_vs_vector(&vs, vids_vector)
    if ret != 0:
        igraph_vector_destroy(&result)
        return -1

    ret = igraph_transitivity_barrat(
        graph,
        &result,
        vs,
        edge_weights,
        IGRAPH_TRANSITIVITY_ZERO
    )

    igraph_vs_destroy(&vs)

    if ret != 0:
        bf_nogil_logf_notime(
            b"LEIDEN",
            "[ERROR] igraph_transitivity_barrat failed (ret=%d) for %u vertices\n",
            ret,
            n_members,
        )
        igraph_vector_destroy(&result)
        return -1

    for i in range(n_members):
        clustering_coefficients[i] = <float>get_vector_element(&result, <igraph_integer_t>i)

    igraph_vector_destroy(&result)
    return 0


# Helper: process a single community (runs nogil)
cdef int process_single_community_global(
    igraph_t* ig_graph,
    igraph_vector_t* ig_weights,
    uint32_t* member_list,
    uint32_t n_members,
    char* keep_flag,
    float* out_comm_cc,
    float* out_individual_ccs,
    float* out_threshold,
    uint32_t* out_node_degrees,
    uint32_t comm_id,
    bint verbose,
    int outlier_method,  # 0 = MAD (default), 1 = IQR, 2 = IFOREST, 3 = LOF, 4 = ZSCORE
    int32_t iforest_n_trees,
    uint32_t iforest_subsample_size,
    double iforest_contamination,
    uint32_t iforest_random_seed,
    uint32_t lof_k,
    double lof_contamination,
    double zscore_threshold,
    ReferenceStats* ref_stats,  # Per-reference statistics for multi-metric outlier detection
    float* out_anomaly_scores,  # Output anomaly scores for each reference (for multi-metric methods)
    uint32_t* exact_connection_counts,  # Real graph metric: exact number of connected neighbors
    double* co_mapping_averages,         # Real graph metric: average co-mapping intensity per read
    uint64_t* max_co_mappings,           # Real graph metric: maximum co-mapping intensity
    double* neighbor_multimap_avg,       # Real graph metric: average multimap rate of neighbors
    float* betweenness_centrality        # Real graph metric: betweenness centrality from igraph
) nogil:
    """
    Compute Barrat clustering coefficients for a community (given as a member_list
    of global vertex ids), apply broken-stick threshold, and set keep_flag for
    kept vertices. Runs nogil and returns 0 on success.

    IMPORTANT: Extracts the induced subgraph for this community to calculate CC
    on the correct topology (only edges within the community).
    """
    cdef int ret
    cdef igraph_vector_int_t tmp_vec
    cdef uint32_t i
    cdef float* comm_clustering = NULL
    cdef float threshold
    cdef uint32_t kept_in_community
    cdef uint32_t best_i
    cdef float max_c
    cdef igraph_t subgraph
    cdef igraph_vector_t sub_weights

    if verbose and not bf_should_log(BF_LOG_LEVEL_TRACE):
        verbose = False

    if n_members == 0:
        return 0

    if n_members == 1:
        keep_flag[member_list[0]] = 1
        return 0

    if n_members == 2:
        keep_flag[member_list[0]] = 1
        keep_flag[member_list[1]] = 1
        return 0

    # Build vertex list for this community
    ret = igraph_vector_int_init(&tmp_vec, <igraph_integer_t>n_members)
    if ret != IGRAPH_SUCCESS:
        return -1

    for i in range(n_members):
        igraph_vector_int_set(&tmp_vec, <igraph_integer_t>i, <igraph_integer_t>member_list[i])

    # Extract induced subgraph for this community only
    # This ensures CC is calculated on the correct topology (intra-community edges only)
    ret = extract_subgraph_with_weights(ig_graph, ig_weights, &tmp_vec, &subgraph, &sub_weights)
    if ret != 0:
        if verbose:
            bf_nogil_logf_notime(
                b"LEIDEN",
                "[WARNING] Failed to extract subgraph for community %u (size=%u), keeping all members\n",
                comm_id,
                n_members,
            )
        igraph_vector_int_destroy(&tmp_vec)
        # Auto-keep all members on error
        for i in range(n_members):
            keep_flag[member_list[i]] = 1
        return 0

    comm_clustering = <float*>malloc(n_members * sizeof(float))
    if not comm_clustering:
        igraph_destroy(&subgraph)
        igraph_vector_destroy(&sub_weights)
        igraph_vector_int_destroy(&tmp_vec)
        return -1

    # Calculate CC on the subgraph (node IDs are now 0 to n_members-1)
    # Build a new vertex selector for the subgraph (all vertices)
    igraph_vector_int_destroy(&tmp_vec)
    ret = igraph_vector_int_init(&tmp_vec, <igraph_integer_t>n_members)
    if ret != IGRAPH_SUCCESS:
        free(comm_clustering)
        igraph_destroy(&subgraph)
        igraph_vector_destroy(&sub_weights)
        return -1

    for i in range(n_members):
        igraph_vector_int_set(&tmp_vec, <igraph_integer_t>i, <igraph_integer_t>i)

    ret = calculate_weighted_clustering_coefficients_for_vids(&subgraph, &sub_weights, &tmp_vec, comm_clustering, n_members)

    if ret != 0:
        # CC calculation failed - cleanup and return
        igraph_vector_int_destroy(&tmp_vec)
        igraph_destroy(&subgraph)
        igraph_vector_destroy(&sub_weights)
        if verbose:
            bf_nogil_logf_notime(
                b"LEIDEN",
                "[WARNING] CC calculation failed for community %u (size=%u), keeping all members\n",
                comm_id,
                n_members,
            )
        for i in range(n_members):
            keep_flag[member_list[i]] = 1
        free(comm_clustering)
        return 0

    # Calculate node degrees in the subgraph (same topology as CC calculation)
    cdef igraph_vector_int_t subgraph_degrees
    cdef igraph_vs_t vs_sub_all
    if out_node_degrees is not NULL:
        ret = igraph_vector_int_init(&subgraph_degrees, 0)
        if ret == IGRAPH_SUCCESS:
            ret = igraph_vs_all(&vs_sub_all)
            if ret == IGRAPH_SUCCESS:
                ret = igraph_degree(&subgraph, &subgraph_degrees, vs_sub_all, IGRAPH_ALL, 0)
                igraph_vs_destroy(&vs_sub_all)
                if ret == IGRAPH_SUCCESS:
                    # Export subgraph degrees (map subgraph IDs back to full graph IDs)
                    for i in range(n_members):
                        out_node_degrees[member_list[i]] = <uint32_t>get_vector_int_element(&subgraph_degrees, <igraph_integer_t>i)
                igraph_vector_int_destroy(&subgraph_degrees)

    # Cleanup subgraph and vertex list
    igraph_vector_int_destroy(&tmp_vec)
    igraph_destroy(&subgraph)
    igraph_vector_destroy(&sub_weights)

    # Calculate statistical outlier threshold for CC values
    # Supports multiple methods: MAD, IQR (univariate), or IFOREST/LOF (multivariate)
    # INTERPRETATION: CC >= threshold → KEEP (cohesive), CC < threshold → REMOVE (hub/contamination)

    cdef OutlierDetectionResult* outlier_result = NULL
    cdef ReferenceFeatures* features = NULL
    cdef uint32_t i_member
    cdef bint* is_outlier_flags = NULL

    # For multi-metric methods (IFOREST, LOF), use the new outlier detection module
    if outlier_method >= 2:  # IFOREST (2), LOF (3), or ZSCORE (4)
        # Allocate feature vectors for community members
        features = <ReferenceFeatures*>malloc(n_members * sizeof(ReferenceFeatures))
        if features == NULL:
            threshold = -1.0e10  # Fallback: keep everything
        else:
            # Extract features from ref_stats for each community member
            for i_member in range(n_members):
                ref_idx = member_list[i_member]
                features[i_member].reference_idx = ref_idx
                features[i_member].clustering_coefficient = <double>comm_clustering[i_member]
                features[i_member].node_degree = <double>out_node_degrees[ref_idx]

                # Extract additional metrics - use REAL graph metrics instead of approximations
                if ref_stats != NULL:
                    features[i_member].multimap_percentage = (
                        <double>ref_stats[ref_idx].repeat_reads / <double>ref_stats[ref_idx].total_reads * 100.0
                        if ref_stats[ref_idx].total_reads > 0 else 0.0
                    )
                else:
                    features[i_member].multimap_percentage = 0.0

                # Use REAL graph metrics (passed as function parameters)
                if exact_connection_counts != NULL:
                    features[i_member].connected_neighbors = <double>exact_connection_counts[ref_idx]
                else:
                    features[i_member].connected_neighbors = 0.0

                if max_co_mappings != NULL:
                    features[i_member].max_comappings = <double>max_co_mappings[ref_idx]
                else:
                    features[i_member].max_comappings = 0.0

                if co_mapping_averages != NULL:
                    features[i_member].avg_comappings_per_read = co_mapping_averages[ref_idx]
                else:
                    features[i_member].avg_comappings_per_read = 0.0

                if neighbor_multimap_avg != NULL:
                    features[i_member].neighbor_multimap_rate = neighbor_multimap_avg[ref_idx] * 100.0
                else:
                    features[i_member].neighbor_multimap_rate = 0.0

                # Use betweenness centrality from igraph
                if betweenness_centrality != NULL:
                    features[i_member].betweenness_centrality = <double>betweenness_centrality[ref_idx]
                else:
                    features[i_member].betweenness_centrality = 0.0

            # Run multi-metric outlier detection
            if outlier_method == OUTLIER_IFOREST:
                outlier_result = detect_outliers_iforest_c(
                    features,
                    n_members,
                    <uint32_t>iforest_n_trees,
                    iforest_subsample_size,
                    iforest_contamination,
                    iforest_random_seed
                )
            elif outlier_method == OUTLIER_LOF:
                outlier_result = detect_outliers_lof_c(
                    features,
                    n_members,
                    lof_k,
                    lof_contamination
                )
            else:
                outlier_result = detect_outliers_multivariate_c(
                    features,
                    n_members,
                    <OutlierMethod>outlier_method,
                    iforest_contamination,
                    iforest_random_seed
                )

            if outlier_result != NULL:
                # Use outlier flags to set keep_flag
                is_outlier_flags = outlier_result.is_outlier
                threshold = <float>outlier_result.threshold

                # Store anomaly scores for all community members (for TSV export)
                if out_anomaly_scores != NULL:
                    for i_member in range(n_members):
                        ref_idx = member_list[i_member]
                        out_anomaly_scores[ref_idx] = <float>outlier_result.anomaly_scores[i_member]
            else:
                # Fallback to MAD if multi-metric detection fails
                threshold = calculate_statistical_outlier_threshold_mad(
                    comm_clustering, n_members,
                    (verbose and bf_should_log(BF_LOG_LEVEL_TRACE)) and n_members > 50
                )

            free(features)
    elif outlier_method == 1:
        # IQR method (univariate on CC only)
        threshold = calculate_statistical_outlier_threshold_iqr(comm_clustering, n_members, (verbose and bf_should_log(BF_LOG_LEVEL_TRACE)) and n_members > 50)
    else:
        # MAD method (default, outlier_method == 0, univariate on CC only)
        threshold = calculate_statistical_outlier_threshold_mad(comm_clustering, n_members, (verbose and bf_should_log(BF_LOG_LEVEL_TRACE)) and n_members > 50)

    # Audit logging: report per-community chosen threshold and community size.
    # printf is safe to call in nogil contexts.
    if verbose and bf_should_log(BF_LOG_LEVEL_TRACE):
        bf_nogil_logf_notime(
            b"LEIDEN",
            "[COMMUNITY] id=%u size=%u chosen_cc_threshold=%.6f\n",
            comm_id,
            n_members,
            threshold,
        )

    # Apply threshold: KEEP references with CC >= threshold (cohesive/informative)
    # REMOVE references with CC < threshold (hub/star pattern, likely contamination)
    # For multi-metric methods, use outlier flags instead of threshold
    kept_in_community = 0
    if is_outlier_flags != NULL:
        # Multi-metric method: use outlier flags (outlier=TRUE means REMOVE, inlier=FALSE means KEEP)
        for i in range(n_members):
            if not is_outlier_flags[i]:  # is_outlier[i] == FALSE → inlier → KEEP
                keep_flag[member_list[i]] = 1
                kept_in_community += 1
    else:
        # Univariate method: use threshold on clustering coefficient
        for i in range(n_members):
            if comm_clustering[i] >= threshold:
                keep_flag[member_list[i]] = 1
                kept_in_community += 1

    if kept_in_community == 0:
        # Safety: ensure at least one reference kept (highest CC = most cohesive)
        best_i = 0
        max_c = comm_clustering[0]
        for i in range(1, n_members):
            if comm_clustering[i] > max_c:
                max_c = comm_clustering[i]
                best_i = i
        keep_flag[member_list[best_i]] = 1
        kept_in_community = 1

    # Audit logging: report how many were kept/removed in this community
    if verbose and bf_should_log(BF_LOG_LEVEL_TRACE):
        bf_nogil_logf_notime(
            b"LEIDEN",
            "[COMMUNITY] id=%u size=%u threshold=%.6f kept=%u removed=%u\n",
            comm_id,
            n_members,
            threshold,
            kept_in_community,
            <uint32_t>(n_members - kept_in_community),
        )

    # compute and export average CC for this community
    cdef float avg_cc = 0.0
    for i in range(n_members):
        avg_cc += comm_clustering[i]
    avg_cc /= <float>n_members
    if out_comm_cc is not NULL:
        out_comm_cc[0] = avg_cc

    # Export individual CC values for each member (using full graph IDs)
    if out_individual_ccs is not NULL:
        for i in range(n_members):
            out_individual_ccs[member_list[i]] = comm_clustering[i]

    # Export threshold value for all members of this community
    if out_threshold is not NULL:
        out_threshold[0] = threshold

    # Cleanup multi-metric outlier detection result
    if outlier_result != NULL:
        free_outlier_result(outlier_result)

    free(comm_clustering)
    return 0


# ==============================================================================
# SUBGRAPH EXTRACTION
# ==============================================================================

cdef int extract_subgraph_with_weights(
    igraph_t* full_graph,
    igraph_vector_t* full_weights,
    igraph_vector_int_t* vertex_ids,
    igraph_t* subgraph_out,
    igraph_vector_t* subgraph_weights_out
) except -1 nogil:
    """
    Extract induced subgraph and corresponding edge weights.
    Returns 0 on success, -1 on failure.

    IMPORTANT: Uses igraph_induced_subgraph to create the subgraph, which
    preserves vertex order and creates a proper induced subgraph (all edges
    between the specified vertices). Then extracts edge weights by iterating
    over the subgraph edges and looking up their weights in the original graph.
    """
    cdef int ret
    cdef igraph_integer_t num_vertices = igraph_vector_int_size(vertex_ids)
    cdef igraph_integer_t i, subgraph_ecount
    cdef igraph_vs_t vs
    cdef igraph_integer_t from_id, to_id, eid
    cdef igraph_integer_t global_vid_from, global_vid_to

    if num_vertices == 0:
        return -1

    # Create vertex selector from vertex_ids
    ret = igraph_vs_vector(&vs, vertex_ids)
    if ret != IGRAPH_SUCCESS:
        return -1

    # Create induced subgraph (includes all edges between specified vertices)
    ret = igraph_induced_subgraph(full_graph, subgraph_out, vs, IGRAPH_SUBGRAPH_AUTO)
    igraph_vs_destroy(&vs)

    if ret != IGRAPH_SUCCESS:
        return -1

    # Get number of edges in subgraph
    subgraph_ecount = igraph_ecount(subgraph_out)

    # Initialize weights vector for subgraph
    ret = igraph_vector_init(subgraph_weights_out, subgraph_ecount)
    if ret != IGRAPH_SUCCESS:
        igraph_destroy(subgraph_out)
        return -1

    # For each edge in the subgraph, look up its weight in the original graph
    # The vertices in the subgraph are remapped (0..num_vertices-1), but we
    # can map them back to the original graph using vertex_ids
    for i in range(subgraph_ecount):
        # Get edge endpoints in subgraph (local vertex IDs: 0..num_vertices-1)
        ret = igraph_edge(subgraph_out, <igraph_integer_t>i, &from_id, &to_id)
        if ret != IGRAPH_SUCCESS:
            igraph_vector_destroy(subgraph_weights_out)
            igraph_destroy(subgraph_out)
            return -1

        # Map subgraph vertex IDs back to original graph vertex IDs
        global_vid_from = get_vector_int_element(vertex_ids, from_id)
        global_vid_to = get_vector_int_element(vertex_ids, to_id)

        # Get edge ID in original graph
        ret = igraph_get_eid(full_graph, &eid, global_vid_from, global_vid_to, IGRAPH_UNDIRECTED, False)
        if ret != IGRAPH_SUCCESS:
            # Edge doesn't exist in original graph (shouldn't happen for induced subgraph)
            igraph_vector_destroy(subgraph_weights_out)
            igraph_destroy(subgraph_out)
            return -1

        # Copy weight from original graph to subgraph weights vector
        (<igraph_real_t**>subgraph_weights_out)[0][i] = get_vector_element(full_weights, eid)

    return 0


# ==============================================================================
# BROKEN STICK MODEL
# ==============================================================================

cdef int _float_compare_ascending(const void* a, const void* b) noexcept nogil:
    """Compare function for qsort (ascending order)."""
    cdef float fa = (<float*>a)[0]
    cdef float fb = (<float*>b)[0]
    if fa < fb:
        return -1
    elif fa > fb:
        return 1
    else:
        return 0


cdef float calculate_statistical_outlier_threshold_mad(float* values, uint32_t n, bint verbose) nogil:
    """
    Calculate threshold for identifying low-CC hubs using Modified Z-Score (MAD method).

    APPROACH:
    Uses Median Absolute Deviation (MAD) for robust outlier detection.
    The modified z-score is: 0.6745 * (x - median) / MAD
    Values with |modified_z_score| > 3.5 are considered outliers.

    FILTERING INTERPRETATION:
    - CC >= threshold → cohesive/well-connected → informative reference → KEEP
    - CC < threshold → hub/star pattern (statistical outlier) → promiscuous/contamination → REMOVE

    High CC means the reference's neighbors are also connected to each other (clique/triangle).
    Low CC means the reference connects otherwise unrelated references (hub/star, likely contaminant).

    ADVANTAGES:
    - Only removes TRUE statistical outliers (may remove 0 references if data is clean!)
    - Robust to extreme outliers (unlike standard deviation)
    - Non-parametric (doesn't assume normal distribution)
    - Well-established in statistics literature

    Args:
        values: Array of clustering coefficients (range [0,1])
        n: Number of values
        verbose: Print diagnostic info

    Returns:
        Threshold value (REMOVE if CC < threshold, KEEP if CC >= threshold)
    """
    cdef uint32_t i
    cdef float* sorted_values = <float*>malloc(n * sizeof(float))
    cdef float* abs_deviations = NULL
    cdef float threshold
    cdef float median_cc, mad, modified_z_threshold
    cdef uint32_t median_idx
    cdef uint32_t num_outliers

    if not sorted_values:
        return 0.0  # No filtering if allocation fails

    # Copy and sort values ascending
    for i in range(n):
        sorted_values[i] = values[i]
    qsort(sorted_values, n, sizeof(float), _float_compare_ascending)

    # Calculate median
    median_idx = n / 2
    if n % 2 == 0 and n > 1:
        median_cc = (sorted_values[median_idx - 1] + sorted_values[median_idx]) / 2.0
    else:
        median_cc = sorted_values[median_idx]

    # Calculate absolute deviations from median
    abs_deviations = <float*>malloc(n * sizeof(float))
    if not abs_deviations:
        free(sorted_values)
        return 0.0

    for i in range(n):
        abs_deviations[i] = sorted_values[i] - median_cc
        if abs_deviations[i] < 0:
            abs_deviations[i] = -abs_deviations[i]

    # Sort absolute deviations to find MAD
    qsort(abs_deviations, n, sizeof(float), _float_compare_ascending)

    # MAD is the median of absolute deviations
    if n % 2 == 0 and n > 1:
        mad = (abs_deviations[median_idx - 1] + abs_deviations[median_idx]) / 2.0
    else:
        mad = abs_deviations[median_idx]

    free(abs_deviations)

    # Calculate key percentiles for diagnostics
    cdef float q25, q75, iqr
    cdef uint32_t q25_idx = n / 4
    cdef uint32_t q75_idx = 3 * n / 4
    q25 = sorted_values[q25_idx]
    q75 = sorted_values[q75_idx]
    iqr = q75 - q25

    if verbose:
        bf_nogil_logf_notime(
            b"LEIDEN",
            "    CC distribution: Q25=%.3f, Median=%.3f, Q75=%.3f, IQR=%.3f, MAD=%.3f\n",
            q25,
            median_cc,
            q75,
            iqr,
            mad,
        )
        bf_nogil_logf_notime(b"LEIDEN", "    Sample size: %u\n", n)

    # Modified Z-score threshold: 3.5 is standard for outlier detection
    # We want to find LOW outliers (hubs), so we look at the lower tail
    # threshold = median - 3.5 * MAD / 0.6745
    if mad > 0.0:
        modified_z_threshold = 3.5
        threshold = median_cc - (modified_z_threshold * mad / 0.6745)

        # Ensure threshold is in valid range [0, 1]
        if threshold < 0.0:
            threshold = 0.0
        if threshold > median_cc:
            threshold = median_cc

        # Count how many values fall below threshold (outliers)
        num_outliers = 0
        for i in range(n):
            if sorted_values[i] < threshold:
                num_outliers += 1

        if verbose:
            bf_nogil_logf_notime(
                b"LEIDEN",
                "    [MAD] Outlier threshold: %.6f (MAD-based, modified_z > 3.5)\n",
                threshold,
            )
            bf_nogil_logf_notime(
                b"LEIDEN",
                "    [MAD] Will REMOVE %u outliers (%.1f%% of community)\n",
                num_outliers,
                100.0 * <float>num_outliers / <float>n,
            )
    else:
        # MAD is 0 (all values identical) - no outliers to remove
        threshold = median_cc
        if verbose:
            bf_nogil_logf_notime(
                b"LEIDEN",
                "    [MAD] MAD=0 (all values identical), no outliers to remove\n",
            )

    free(sorted_values)
    return threshold


cdef float calculate_statistical_outlier_threshold_iqr(float* values, uint32_t n, bint verbose) nogil:
    """
    Calculate threshold for identifying low-CC hubs using IQR method.

    APPROACH:
    Uses Interquartile Range (IQR) for outlier detection.
    Lower bound = Q1 - 1.5 * IQR
    Upper bound = Q3 + 1.5 * IQR (not used here, we only care about low outliers)

    FILTERING INTERPRETATION:
    - CC >= threshold → cohesive/well-connected → informative reference → KEEP
    - CC < threshold → hub/star pattern (statistical outlier) → promiscuous/contamination → REMOVE

    ADVANTAGES:
    - Only removes TRUE statistical outliers (may remove 0 references if data is clean!)
    - Standard boxplot outlier method
    - Non-parametric (doesn't assume normal distribution)
    - Simple and interpretable

    Args:
        values: Array of clustering coefficients (range [0,1])
        n: Number of values
        verbose: Print diagnostic info

    Returns:
        Threshold value (REMOVE if CC < threshold, KEEP if CC >= threshold)
    """
    cdef uint32_t i
    cdef float* sorted_values = <float*>malloc(n * sizeof(float))
    cdef float threshold
    cdef float q25, q75, iqr, lower_bound
    cdef uint32_t q25_idx, q75_idx, median_idx
    cdef float median_cc
    cdef uint32_t num_outliers

    if not sorted_values:
        return 0.0  # No filtering if allocation fails

    # Copy and sort values ascending
    for i in range(n):
        sorted_values[i] = values[i]
    qsort(sorted_values, n, sizeof(float), _float_compare_ascending)

    # Calculate quartiles
    q25_idx = n / 4
    q75_idx = 3 * n / 4
    median_idx = n / 2

    q25 = sorted_values[q25_idx]
    q75 = sorted_values[q75_idx]
    median_cc = sorted_values[median_idx]
    iqr = q75 - q25

    # Lower outlier bound: Q1 - 1.5 * IQR
    lower_bound = q25 - 1.5 * iqr

    # Ensure threshold is in valid range [0, 1]
    threshold = lower_bound
    if threshold < 0.0:
        threshold = 0.0
    if threshold > 1.0:
        threshold = 1.0

    # Count how many values fall below threshold (outliers)
    num_outliers = 0
    for i in range(n):
        if sorted_values[i] < threshold:
            num_outliers += 1

    if verbose:
        bf_nogil_logf_notime(
            b"LEIDEN",
            "    CC distribution: Q25=%.3f, Median=%.3f, Q75=%.3f, IQR=%.3f\n",
            q25,
            median_cc,
            q75,
            iqr,
        )
        bf_nogil_logf_notime(b"LEIDEN", "    Sample size: %u\n", n)
        bf_nogil_logf_notime(
            b"LEIDEN",
            "    [IQR] Outlier threshold: %.6f (Q1 - 1.5*IQR = %.3f - 1.5*%.3f)\n",
            threshold,
            q25,
            iqr,
        )
        bf_nogil_logf_notime(
            b"LEIDEN",
            "    [IQR] Will REMOVE %u outliers (%.1f%% of community)\n",
            num_outliers,
            100.0 * <float>num_outliers / <float>n,
        )

    free(sorted_values)
    return threshold


cdef float calculate_broken_stick_threshold(float* values, uint32_t n, bint verbose) nogil:
    """
    Statistical outlier detection using MAD (Median Absolute Deviation).

    Uses Modified Z-Score with MAD for robust outlier detection.
    Only removes TRUE statistical outliers (may remove 0 references if data is clean).

    This function maintains the old name for backward compatibility but uses
    the new MAD-based statistical approach instead of percentiles.
    """
    return calculate_statistical_outlier_threshold_mad(values, n, verbose)


# ==============================================================================
# COMPONENT INFO STRUCTURE
# ==============================================================================

cdef struct ComponentInfo:
    uint32_t component_id
    uint32_t* members
    uint32_t n_members
    uint32_t refs_removed
    double time_total
    double time_extract
    double time_barrat
    double time_leiden
    double time_communities


# ==============================================================================
# LEIDEN RESULTS STRUCTURE
# ==============================================================================

cdef struct LeidenResults:
    uint32_t num_nodes  # number of nodes processed
    char* keep_flag  # 1 = keep, 0 = remove
    uint32_t* community_membership  # community ID for each reference
    float* community_cc_values  # Average CC value for the community each reference belongs to
    float* individual_cc_values  # Individual CC value for each reference (Barrat's method)
    float* cc_threshold_values  # Broken-stick threshold used for each reference's community
    float* anomaly_scores  # Anomaly score for each reference (for multi-metric methods like Isolation Forest)
    uint32_t* node_degree  # Node degree (number of edges) - already in ReferencePattern, kept here for convenience


# ==============================================================================
# MAIN LEIDEN CLUSTERING ENTRY POINT
# ==============================================================================

cdef LeidenResults* leiden_clustering(WeightedGraph* graph,
                             ReferenceStats* ref_stats,
                             double resolution,
                             int32_t max_iterations,
                             bint verbose,
                             int32_t thread_count,
                             int outlier_method,
                             int32_t iforest_n_trees,
                             uint32_t iforest_subsample_size,
                             double iforest_contamination,
                             uint32_t iforest_random_seed,
                             uint32_t lof_k,
                             double lof_contamination,
                             double zscore_threshold,
                             uint32_t* exact_connection_counts,
                             double* co_mapping_averages,
                             uint64_t* max_co_mappings,
                             double* neighbor_multimap_avg,
                             uint32_t array_size) except NULL nogil:
    """
    Run Leiden clustering with clustering coefficient filtering.
    
    Returns LeidenResults struct containing:
    - keep_flag: 1 = keep, 0 = remove
    - community_membership: community ID for each reference
    - community_cc_values: CC value for the community each reference belongs to
    
    Workflow:
    1. Get cached igraph from WeightedGraph
    2. Find connected components
    3. Auto-keep all singletons (no processing needed)
    4. For each multi-node component, run Leiden + CC filtering (PARALLEL)
    5. Return LeidenResults struct
    """
    cdef igraph_t* ig_graph
    cdef igraph_vector_t* ig_weights
    cdef igraph_vector_int_t component_membership, component_sizes
    cdef igraph_integer_t num_components
    cdef char* keep_flag = NULL
    cdef ComponentInfo* components = NULL
    cdef uint32_t** component_members_lists = NULL
    cdef uint32_t* component_counts = NULL
    cdef uint32_t i, k, ref_idx, comm_id
    cdef int ret
    cdef uint32_t idx
    cdef int effective_threads
    cdef uint32_t singletons_count = 0
    cdef uint32_t multi_node_count = 0
    cdef uint32_t* multi_node_indices = NULL
    cdef uint32_t multi_idx
    # Stats variables for multi-node component sizes
    cdef uint32_t min_size, max_size, sidx, cur_size
    cdef double sum_sizes
    # Timing aggregation temporaries (declare at function top for Cython)
    cdef double total_time, total_extract, total_barrat, total_leiden, total_communities
    cdef double t_start, t_end
    cdef uint32_t ci
    cdef uint32_t key, key_size
    cdef int j
    # Additional C-level locals hoisted here to satisfy Cython's requirement
    # that all cdef declarations appear before executable statements.
    cdef igraph_vector_int_t global_membership
    cdef igraph_integer_t num_communities
    cdef igraph_vector_int_t community_sizes
    cdef igraph_vector_int_t community_membership_tmp
    cdef igraph_integer_t max_m
    cdef igraph_integer_t tmp_m
    cdef uint32_t num_comms_u
    cdef uint32_t comm_idx
    cdef uint32_t n_members_in_comm
    cdef uint32_t actual_num_comms
    cdef uint32_t* community_member_counts
    cdef uint32_t** community_member_lists
    cdef uint32_t member_pos
    cdef uint32_t large_comm_count
    cdef uint32_t singleton_comm_count
    cdef uint32_t pair_comm_count
    cdef uint32_t* active_communities
    cdef uint32_t ac_pos
    cdef uint32_t comms_to_process
    cdef igraph_vector_int_t vids_vec
    cdef igraph_vector_int_t tmp_vec
    cdef double t_leiden_total
    cdef igraph_integer_t v, m
    cdef igraph_integer_t gm

    if verbose and not bf_should_log(BF_LOG_LEVEL_TRACE):
        verbose = False
    cdef igraph_t subgraph
    cdef igraph_vector_t sub_weights
    cdef igraph_vector_int_t vids_to_leiden
    cdef igraph_vector_int_t membership_sel
    cdef igraph_vs_t vs_sel
    cdef igraph_integer_t n_selected
    cdef float* per_comm_avg
    cdef uint32_t n_members
    cdef float comm_cc
    cdef uint32_t total_assigned
    
    if verbose:
        bf_nogil_logf_notime(b"LEIDEN", "\n=== LEIDEN CLUSTERING WITH CC FILTERING (PARALLEL) ===\n")
        bf_nogil_logf_notime(b"LEIDEN", "Resolution: %.3f\n", resolution)
        bf_nogil_logf_notime(b"LEIDEN", "Max iterations: %d\n", max_iterations)
        bf_nogil_logf_notime(b"LEIDEN", "Thread count: %d\n", thread_count)
    
    # Get cached igraph and weights pointers
    if not graph.igraph_handle:
        bf_nogil_logf_notime(b"LEIDEN", "ERROR: No cached igraph found in WeightedGraph.igraph_handle!\n")
        return NULL
    
    if not graph.weights_handle:
        bf_nogil_logf_notime(b"LEIDEN", "ERROR: No cached weights found in WeightedGraph.weights_handle!\n")
        return NULL
    
    ig_graph = <igraph_t*>graph.igraph_handle
    ig_weights = <igraph_vector_t*>graph.weights_handle
    
    if verbose:
        bf_nogil_logf_notime(
            b"LEIDEN",
            "Graph: %ld vertices, %ld edges\n",
            <long>igraph_vcount(ig_graph),
            <long>igraph_ecount(ig_graph),
        )
        bf_nogil_logf_notime(b"LEIDEN", "Finding connected components...\n")
    
    # Find connected components
    ret = igraph_vector_int_init(&component_membership, graph.num_nodes)
    if ret != IGRAPH_SUCCESS:
        bf_nogil_logf_notime(b"LEIDEN", "ERROR: Failed to initialize component membership vector\n")
        return NULL
    
    ret = igraph_vector_int_init(&component_sizes, 0)
    if ret != IGRAPH_SUCCESS:
        igraph_vector_int_destroy(&component_membership)
        bf_nogil_logf_notime(b"LEIDEN", "ERROR: Failed to initialize component sizes vector\n")
        return NULL
    
    ret = igraph_clusters(ig_graph, &component_membership, &component_sizes,
                         &num_components, IGRAPH_WEAK)
    if ret != IGRAPH_SUCCESS:
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"LEIDEN", "ERROR: Failed to find connected components\n")
        return NULL
    
    if verbose:
        bf_nogil_logf_notime(b"LEIDEN", "Found %ld connected components\n", <long>num_components)
    
    # Allocate LeidenResults struct
    cdef LeidenResults* results = <LeidenResults*>malloc(sizeof(LeidenResults))
    if not results:
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"LEIDEN", "ERROR: Failed to allocate LeidenResults\n")
        return NULL
    
    results.num_nodes = graph.num_nodes
    
    # Allocate keep_flag
    results.keep_flag = <char*>calloc(graph.num_nodes, sizeof(char))
    if not results.keep_flag:
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"LEIDEN", "ERROR: Failed to allocate keep_flag\n")
        return NULL
    
    # Allocate community membership array
    results.community_membership = <uint32_t*>malloc(graph.num_nodes * sizeof(uint32_t))
    if not results.community_membership:
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"LEIDEN", "ERROR: Failed to allocate community_membership\n")
        return NULL

    # Allocate component membership array
    results.component_membership = <uint32_t*>malloc(graph.num_nodes * sizeof(uint32_t))
    if not results.component_membership:
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"LEIDEN", "ERROR: Failed to allocate component_membership\n")
        return NULL

    # Allocate node degree array
    results.node_degree = <uint32_t*>malloc(graph.num_nodes * sizeof(uint32_t))
    if not results.node_degree:
        free(results.node_degree)
        free(results.component_membership)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"LEIDEN", "ERROR: Failed to allocate node_degree\n")
        return NULL

    # Allocate community CC values array
    results.community_cc_values = <float*>malloc(graph.num_nodes * sizeof(float))
    if not results.community_cc_values:
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"LEIDEN", "ERROR: Failed to allocate community_cc_values\n")
        return NULL

    # Allocate individual CC values array
    results.individual_cc_values = <float*>malloc(graph.num_nodes * sizeof(float))
    if not results.individual_cc_values:
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"LEIDEN", "ERROR: Failed to allocate individual_cc_values\n")
        return NULL

    # Allocate CC threshold values array
    results.cc_threshold_values = <float*>malloc(graph.num_nodes * sizeof(float))
    if not results.cc_threshold_values:
        free(results.individual_cc_values)
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"LEIDEN", "ERROR: Failed to allocate cc_threshold_values\n")
        return NULL

    # Allocate anomaly scores array (for multi-metric methods like Isolation Forest)
    results.anomaly_scores = <float*>calloc(graph.num_nodes, sizeof(float))
    if not results.anomaly_scores:
        free(results.cc_threshold_values)
        free(results.individual_cc_values)
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"LEIDEN", "ERROR: Failed to allocate anomaly_scores\n")
        return NULL

    # Allocate betweenness centrality array
    results.betweenness_centrality = <float*>calloc(graph.num_nodes, sizeof(float))
    if not results.betweenness_centrality:
        free(results.anomaly_scores)
        free(results.cc_threshold_values)
        free(results.individual_cc_values)
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"LEIDEN", "ERROR: Failed to allocate anomaly_scores\n")
        return NULL

    # Copy component membership from igraph vector to results
    for ref_idx in range(graph.num_nodes):
        results.component_membership[ref_idx] = <uint32_t>get_vector_int_element(&component_membership, ref_idx)

    # Initialize node degrees to 0 (will be filled during community processing)
    # Node degrees will be calculated from community subgraphs, not the full graph
    for ref_idx in range(graph.num_nodes):
        results.node_degree[ref_idx] = 0

    # Calculate betweenness centrality for all nodes using igraph
    # This is a global graph metric (not community-specific)
    if verbose:
        bf_nogil_logf_notime(b"LEIDEN", "Calculating betweenness centrality...\n")

    cdef igraph_vector_t betweenness_vec
    cdef igraph_vs_t vs_all
    ret = igraph_vector_init(&betweenness_vec, 0)
    if ret == IGRAPH_SUCCESS:
        ret = igraph_vs_all(&vs_all)
        if ret == IGRAPH_SUCCESS:
            # Calculate betweenness centrality (directed=False, weights=ig_weights)
            ret = igraph_betweenness(ig_graph, &betweenness_vec, vs_all, 0, ig_weights)
            igraph_vs_destroy(&vs_all)

            if ret == IGRAPH_SUCCESS:
                # Copy betweenness values to results
                for ref_idx in range(graph.num_nodes):
                    results.betweenness_centrality[ref_idx] = <float>VECTOR(betweenness_vec)[ref_idx]

                if verbose:
                    bf_nogil_logf_notime(b"LEIDEN", "  Betweenness centrality calculated for %u nodes\n", graph.num_nodes)
            else:
                if verbose:
                    bf_nogil_logf_notime(b"LEIDEN", "[WARNING] Betweenness calculation failed, using zeros\n")

            igraph_vector_destroy(&betweenness_vec)
        else:
            if verbose:
                bf_nogil_logf_notime(b"LEIDEN", "[WARNING] Failed to create vertex selector for betweenness\n")
    else:
        if verbose:
            bf_nogil_logf_notime(b"LEIDEN", "[WARNING] Failed to initialize betweenness vector\n")

    # Count members per component and identify singletons
    component_counts = <uint32_t*>calloc(num_components, sizeof(uint32_t))
    if not component_counts:
        free(results.node_degree)
        free(results.component_membership)
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        return NULL

    for ref_idx in range(graph.num_nodes):
        comm_id = get_vector_int_element(&component_membership, ref_idx)
        if comm_id < num_components:
            component_counts[comm_id] += 1
    
    # Separate singletons from multi-node components
    for k in range(num_components):
        if component_counts[k] == 1:
            singletons_count += 1
        elif component_counts[k] > 1:
            multi_node_count += 1
    
    if verbose:
        bf_nogil_logf_notime(b"LEIDEN", "  Singletons: %u (auto-keep, no processing)\n", singletons_count)
        bf_nogil_logf_notime(b"LEIDEN", "  Multi-node components: %u (will process)\n", multi_node_count)
    
    # Auto-keep all singleton nodes
    for ref_idx in range(graph.num_nodes):
        comm_id = get_vector_int_element(&component_membership, ref_idx)
        if comm_id < num_components and component_counts[comm_id] == 1:
            results.keep_flag[ref_idx] = 1

    # Also auto-keep pairs (component size == 2) if we exclude pairs from Leiden
    for ref_idx in range(graph.num_nodes):
        comm_id = get_vector_int_element(&component_membership, ref_idx)
        if comm_id < num_components and component_counts[comm_id] == 2:
            results.keep_flag[ref_idx] = 1
    
    # If no multi-node components, we're done
    if multi_node_count == 0:
        if verbose:
            bf_nogil_logf_notime(b"LEIDEN", "No multi-node components to process!\n")
        free(component_counts)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        return results
    
    # Build list of multi-node component indices
    multi_node_indices = <uint32_t*>malloc(multi_node_count * sizeof(uint32_t))
    if not multi_node_indices:
        free(component_counts)
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        return NULL

    multi_idx = 0
    for k in range(num_components):
        if component_counts[k] > 1:
            multi_node_indices[multi_idx] = k
            multi_idx += 1

    # Compute and print basic size statistics for multi-node components
    # (min, max, mean) to help diagnostics when verbose is enabled.
    if verbose and multi_node_count > 0:
        # Use variables declared at function top (Cython requires cdef at top-level)
        sum_sizes = 0.0

        min_size = component_counts[multi_node_indices[0]]
        max_size = min_size

        for sidx in range(multi_node_count):
            cur_size = component_counts[multi_node_indices[sidx]]
            if cur_size < min_size:
                min_size = cur_size
            if cur_size > max_size:
                max_size = cur_size
            sum_sizes += <double>cur_size

        bf_nogil_logf_notime(
            b"LEIDEN",
            "  Multi-node component sizes: min=%u, max=%u, mean=%.2f\n",
            min_size,
            max_size,
            sum_sizes / <double>multi_node_count,
        )
    
    # New approach: run Leiden once on the full graph, then process communities in parallel
    # Build membership vector for full graph using igraph_community_leiden

    # Allocate global membership vector (size = number of nodes) and set to -1
    # We'll run Leiden only on the non-singleton vertices (faster) and map
    # the resulting membership back into this full-size vector so downstream
    # code can remain unchanged.
    ret = igraph_vector_int_init(&global_membership, graph.num_nodes)
    if ret != IGRAPH_SUCCESS:
        free(component_counts)
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        return NULL

    # initialize all entries to -1 (meaning: not assigned / singleton)
    for i in range(graph.num_nodes):
        igraph_vector_int_set(&global_membership, <igraph_integer_t>i, <igraph_integer_t>-1)

    # Build vertex list excluding singletons (component_counts[comm_id] > 1)

    ret = igraph_vector_int_init(&vids_to_leiden, 0)
    if ret != IGRAPH_SUCCESS:
        igraph_vector_int_destroy(&global_membership)
        free(component_counts)
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        return NULL

    # Select vertices for Leiden: exclude singletons and pairs (component size <= 2)
    for i in range(graph.num_nodes):
        comm_id = get_vector_int_element(&component_membership, i)
        if comm_id >= 0 and comm_id < num_components and component_counts[comm_id] > 2:
            igraph_vector_int_push_back(&vids_to_leiden, <igraph_integer_t>i)

    n_selected = igraph_vector_int_size(&vids_to_leiden)

    if n_selected == 0:
        # Nothing to run Leiden on (all singletons); set num_communities=0
        num_communities = 0
        if verbose:
            bf_nogil_logf_notime(b"LEIDEN", "No non-singleton vertices - skipping Leiden\n")
        igraph_vector_int_destroy(&vids_to_leiden)
    else:
        # Prepare membership vector for selected vertices
        ret = igraph_vector_int_init(&membership_sel, n_selected)
        if ret != IGRAPH_SUCCESS:
            igraph_vector_int_destroy(&vids_to_leiden)
            igraph_vector_int_destroy(&global_membership)
            free(component_counts)
            free(results.community_cc_values)
            free(results.community_membership)
            free(results.keep_flag)
            free(results)
            igraph_vector_int_destroy(&component_membership)
            igraph_vector_int_destroy(&component_sizes)
            return NULL

        ret = extract_subgraph_with_weights(ig_graph, ig_weights, &vids_to_leiden,
                                            &subgraph, &sub_weights)
        if ret != 0:
            if verbose:
                bf_nogil_logf_notime(b"LEIDEN", "ERROR: Failed to extract subgraph for Leiden subset\n")
            igraph_vector_int_destroy(&membership_sel)
            igraph_vector_int_destroy(&vids_to_leiden)
            igraph_vector_int_destroy(&global_membership)
            free(component_counts)
            free(results.community_cc_values)
            free(results.community_membership)
            free(results.keep_flag)
            free(results)
            igraph_vector_int_destroy(&component_membership)
            igraph_vector_int_destroy(&component_sizes)
            return NULL

        # time Leiden on the extracted subgraph
        t_start = clock()
        ret = igraph_community_leiden(&subgraph, &sub_weights, NULL,
                                      <igraph_real_t>resolution, 0.01, False,
                                      <igraph_integer_t>max_iterations,
                                      &membership_sel, &num_communities, NULL)
        t_end = clock()
        t_leiden_total = (t_end - t_start) / <double>CLOCKS_PER_SEC

        # cleanup subgraph resources
        igraph_destroy(&subgraph)
        igraph_vector_destroy(&sub_weights)

        if ret != IGRAPH_SUCCESS:
            if verbose:
                bf_nogil_logf_notime(b"LEIDEN", "ERROR: Leiden failed on extracted subgraph\n")
            igraph_vector_int_destroy(&membership_sel)
            igraph_vector_int_destroy(&vids_to_leiden)
            igraph_vector_int_destroy(&global_membership)
            free(component_counts)
            free(results.community_cc_values)
            free(results.community_membership)
            free(results.keep_flag)
            free(results)
            igraph_vector_int_destroy(&component_membership)
            igraph_vector_int_destroy(&component_sizes)
            return NULL

        if verbose:
            bf_nogil_logf_notime(b"LEIDEN", "Leiden (subset) produced %ld communities in %.3f s\n", <long>num_communities, t_leiden_total)

        # Map membership_sel (length n_selected) back to the full-size global_membership
        for i in range(<igraph_integer_t>n_selected):
            v = get_vector_int_element(&vids_to_leiden, i)
            m = get_vector_int_element(&membership_sel, i)
            igraph_vector_int_set(&global_membership, v, m)

        # Defensive check: recompute actual number of communities from membership_sel
        # Some igraph builds have returned unexpected num_communities values in the
        # wild (larger than the number of vertices selected). Recompute the
        # actual max membership id and adjust num_communities accordingly so
        # subsequent allocations use a safe, correct size.
        max_m = -1
        for i in range(<igraph_integer_t>n_selected):
            tmp_m = get_vector_int_element(&membership_sel, i)
            if tmp_m > max_m:
                max_m = tmp_m

        actual_num_comms = 0
        if max_m >= 0:
            actual_num_comms = <uint32_t>(max_m + 1)
        else:
            actual_num_comms = 0

        if <uint32_t>num_communities != actual_num_comms:
            if verbose:
                bf_nogil_logf_notime(
                    b"LEIDEN",
                    "[LEIDEN] Adjusting num_communities: reported=%ld recomputed=%u\n",
                    <long>num_communities,
                    actual_num_comms,
                )
            # Use the recomputed value for allocations and downstream logic
            num_communities = <igraph_integer_t>actual_num_comms

        igraph_vector_int_destroy(&membership_sel)
        igraph_vector_int_destroy(&vids_to_leiden)

    # Prepare lists for each community (only for communities with size > 0)
    # First count community sizes
    community_member_counts = <uint32_t*>calloc(<size_t>num_communities, sizeof(uint32_t))
    if not community_member_counts:
        igraph_vector_int_destroy(&global_membership)
        free(component_counts)
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        return NULL

    for ref_idx in range(graph.num_nodes):
        comm_id = get_vector_int_element(&global_membership, ref_idx)
        if comm_id >= 0 and comm_id < num_communities:
            community_member_counts[comm_id] += 1

    # Diagnostic: report how many vertices were selected for Leiden
    if verbose:
        bf_nogil_logf_notime(b"LEIDEN", "[LEIDEN] Vertices selected for Leiden: %ld\n", <long>n_selected)

    # Allocate member lists
    community_member_lists = <uint32_t**>malloc(num_communities * sizeof(uint32_t*))
    if not community_member_lists:
        free(community_member_counts)
        igraph_vector_int_destroy(&global_membership)
        free(component_counts)
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        return NULL

    for comm_idx in range(num_communities):
        if community_member_counts[comm_idx] > 0:
            community_member_lists[comm_idx] = <uint32_t*>malloc(community_member_counts[comm_idx] * sizeof(uint32_t))
            if not community_member_lists[comm_idx]:
                # free previous
                for i in range(comm_idx):
                    if community_member_lists[i]: free(community_member_lists[i])
                free(community_member_lists)
                free(community_member_counts)
                igraph_vector_int_destroy(&global_membership)
                free(component_counts)
                free(results.community_cc_values)
                free(results.community_membership)
                free(results.keep_flag)
                free(results)
                igraph_vector_int_destroy(&component_membership)
                igraph_vector_int_destroy(&component_sizes)
                return NULL
        else:
            community_member_lists[comm_idx] = NULL

    # Fill member lists
    # reset counts to use as cursor
    for comm_idx in range(num_communities):
        community_member_counts[comm_idx] = 0

    for ref_idx in range(graph.num_nodes):
        comm_id = get_vector_int_element(&global_membership, ref_idx)
        if comm_id >= 0 and comm_id < num_communities:
            member_pos = community_member_counts[comm_id]
            community_member_lists[comm_id][member_pos] = ref_idx
            community_member_counts[comm_id] += 1
    
    # Initialize CC-related arrays to 0.0
    for ref_idx in range(graph.num_nodes):
        results.community_cc_values[ref_idx] = 0.0
        results.individual_cc_values[ref_idx] = 0.0
        results.cc_threshold_values[ref_idx] = 0.0
    # Build a compact list of communities that need full processing (size > 2)
    large_comm_count = 0
    singleton_comm_count = 0
    pair_comm_count = 0
    for comm_idx in range(num_communities):
        if community_member_counts[comm_idx] == 0:
            continue
        elif community_member_counts[comm_idx] == 1:
            singleton_comm_count += 1
        elif community_member_counts[comm_idx] == 2:
            pair_comm_count += 1
        else:
            large_comm_count += 1

    if verbose:
        bf_nogil_logf_notime(
            b"LEIDEN",
            "Community sizes: singletons=%u, pairs=%u, larger=%u (total communities=%ld)\n",
            singleton_comm_count,
            pair_comm_count,
            large_comm_count,
            <long>num_communities,
        )

    # Diagnostic: sum of all community_member_counts (should equal number of
    # nodes assigned to communities by Leiden; singletons/pairs may not be
    # included in the Leiden selection, so this is a sanity check).
    if verbose:
        total_assigned = 0
        for comm_idx in range(num_communities):
            total_assigned += community_member_counts[comm_idx]
        bf_nogil_logf_notime(b"LEIDEN", "[LEIDEN] Total assigned to communities: %u\n", total_assigned)

    # Allocate per-community buffers (initialized to 0)
    cdef float* per_comm_thresholds = NULL
    if num_communities > 0:
        per_comm_avg = <float*>calloc(num_communities, sizeof(float))
        per_comm_thresholds = <float*>calloc(num_communities, sizeof(float))

    # Process ALL communities in parallel (process_single_community_global handles small ones efficiently)
    effective_threads = thread_count if thread_count > 0 else 1
    if effective_threads > <int>num_communities:
        effective_threads = <int>num_communities

    if verbose:
        bf_nogil_logf_notime(
            b"LEIDEN",
            "Processing %u communities with %d threads...\n",
            <uint32_t>num_communities,
            effective_threads,
        )

    # Process all communities (sizes 1, 2, and >2)
    # The process_single_community_global function handles all sizes efficiently
    for comm_idx in prange(num_communities, nogil=True, schedule='dynamic', num_threads=effective_threads):
        if community_member_counts[comm_idx] == 0:
            continue
        process_single_community_global(ig_graph, ig_weights,
                                        community_member_lists[comm_idx],
                                        community_member_counts[comm_idx],
                                        results.keep_flag,
                                        &per_comm_avg[comm_idx],
                                        results.individual_cc_values,
                                        &per_comm_thresholds[comm_idx],
                                        results.node_degree,
                                        comm_idx,
                                        verbose,
                                        outlier_method,
                                        iforest_n_trees,
                                        iforest_subsample_size,
                                        iforest_contamination,
                                        iforest_random_seed,
                                        lof_k,
                                        lof_contamination,
                                        zscore_threshold,
                                        ref_stats,
                                        results.anomaly_scores,
                                        exact_connection_counts,
                                        co_mapping_averages,
                                        max_co_mappings,
                                        neighbor_multimap_avg,
                                        results.betweenness_centrality)

    # Populate results.community_cc_values and cc_threshold_values from per_comm buffers
    for comm_idx in range(num_communities):
        if community_member_counts[comm_idx] == 0:
            continue
        n_members = community_member_counts[comm_idx]
        if n_members <= 2:
            comm_cc = 0.0
        else:
            comm_cc = per_comm_avg[comm_idx]
        for i in range(n_members):
            ref_idx = community_member_lists[comm_idx][i]
            results.community_cc_values[ref_idx] = comm_cc
            # Store the threshold used for this community
            results.cc_threshold_values[ref_idx] = per_comm_thresholds[comm_idx]

    # Map global membership vector into results.community_membership for all nodes.
    # Nodes that were not part of the Leiden run (singletons) will have -1 in
    # global_membership; represent those as UINT32_MAX (0xFFFFFFFF) in the
    # public results so downstream TSV/BAM writers can encode 'no community'.
    for ref_idx in range(graph.num_nodes):
        gm = get_vector_int_element(&global_membership, <igraph_integer_t>ref_idx)
        if gm >= 0:
            results.community_membership[ref_idx] = <uint32_t>gm
        else:
            results.community_membership[ref_idx] = <uint32_t>0xFFFFFFFF

    # Cleanup community lists
    for comm_idx in range(num_communities):
        if community_member_lists[comm_idx]:
            free(community_member_lists[comm_idx])
    free(community_member_lists)
    free(community_member_counts)
    if per_comm_avg:
        free(per_comm_avg)
    if per_comm_thresholds:
        free(per_comm_thresholds)
    free(component_counts)
    igraph_vector_int_destroy(&component_membership)
    igraph_vector_int_destroy(&component_sizes)
    
    if verbose:
        bf_nogil_logf_notime(b"LEIDEN", "Leiden clustering with CC filtering completed\n")
        bf_nogil_logf_notime(b"LEIDEN", "==========================================================\n\n")
    
    
    return results
