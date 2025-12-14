# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: initializedcheck=False
"""
Network QC metrics for evaluating EM quality based on graph structure.
"""

from libc.math cimport log, log2, sqrt, fabs, fmax, fmin, pow as libc_pow
from libc.stdlib cimport malloc, free, qsort
from libc.string cimport memset, memcpy
from libc.stdint cimport uint32_t, int32_t, uint64_t, int64_t

from bam_filter.processor_network_qc cimport (
    NetworkQCDataset, NetworkQCReference, NetworkQCConfig,
    compute_tax_ambiguity_flags,
)


# ============================================================================
# Constants
# ============================================================================

cdef double LOG2_E = 1.4426950408889634  # 1/ln(2) for converting ln to log2
cdef double EPSILON = 1e-15


# ============================================================================
# Utility: Gini coefficient
# ============================================================================

cdef int _compare_double(const void* a, const void* b) noexcept nogil:
    """Comparison function for qsort on doubles."""
    cdef double da = (<double*>a)[0]
    cdef double db = (<double*>b)[0]
    if da < db:
        return -1
    elif da > db:
        return 1
    return 0


cdef double compute_gini(double* values, uint32_t n) noexcept nogil:
    """
    Compute Gini coefficient for inequality measurement.

    G = (2 * Σ i*x_i) / (n * Σ x_i) - (n+1)/n

    Returns value in [0, 1] where 0 = perfect equality, 1 = maximum inequality.
    """
    if n == 0:
        return 0.0

    # Copy and sort values
    cdef double* sorted_vals = <double*>malloc(n * sizeof(double))
    if sorted_vals == NULL:
        return 0.0

    memcpy(sorted_vals, values, n * sizeof(double))
    qsort(sorted_vals, n, sizeof(double), _compare_double)

    cdef double total = 0.0
    cdef double weighted_sum = 0.0
    cdef uint32_t i

    for i in range(n):
        total += sorted_vals[i]
        weighted_sum += (i + 1) * sorted_vals[i]

    free(sorted_vals)

    if total < EPSILON:
        return 0.0

    return (2.0 * weighted_sum) / (n * total) - (n + 1.0) / n


# ============================================================================
# Utility: Shannon entropy
# ============================================================================

cdef double compute_entropy(double* values, uint32_t n) noexcept nogil:
    """
    Compute Shannon entropy in bits.

    H = -Σ p_i * log2(p_i)

    Input values are normalized internally.
    """
    if n == 0:
        return 0.0

    cdef double total = 0.0
    cdef uint32_t i

    for i in range(n):
        if values[i] > 0:
            total += values[i]

    if total < EPSILON:
        return 0.0

    cdef double entropy = 0.0
    cdef double p

    for i in range(n):
        if values[i] > EPSILON:
            p = values[i] / total
            entropy -= p * log2(p)

    return entropy


# ============================================================================
# Utility: Effective number (inverse Simpson index)
# ============================================================================

cdef double compute_n_eff(double* values, uint32_t n) noexcept nogil:
    """
    Compute effective number of components: N_eff = 1 / Σ p_j²

    This is the inverse Simpson index, indicating diversity.
    """
    if n == 0:
        return 0.0

    cdef double total = 0.0
    cdef double sum_sq = 0.0
    cdef uint32_t i

    for i in range(n):
        if values[i] > 0:
            total += values[i]

    if total < EPSILON:
        return 0.0

    for i in range(n):
        if values[i] > 0:
            sum_sq += (values[i] / total) ** 2

    if sum_sq < EPSILON:
        return <double>n

    return 1.0 / sum_sq


# ============================================================================
# Modularity computation
# ============================================================================

cdef double compute_modularity(
    uint32_t* edge_src,
    uint32_t* edge_dst,
    double* edge_weights,
    uint64_t n_edges,
    int32_t* partition,
    uint32_t n_nodes,
) noexcept nogil:
    """
    Compute modularity Q for a given partition.

    Q = (1/2m) Σ_{ij} [A_{ij} - k_i*k_j/(2m)] * δ(c_i, c_j)

    where:
    - A_{ij} is edge weight between i and j
    - k_i is strength (sum of edge weights) of node i
    - m is total edge weight
    - c_i is community of node i
    """
    if n_edges == 0 or n_nodes == 0:
        return 0.0

    # Compute node strengths and total edge weight
    cdef double* strengths = <double*>malloc(n_nodes * sizeof(double))
    if strengths == NULL:
        return 0.0

    memset(strengths, 0, n_nodes * sizeof(double))

    cdef double total_weight = 0.0
    cdef uint64_t e
    cdef uint32_t src, dst
    cdef double w

    for e in range(n_edges):
        src = edge_src[e]
        dst = edge_dst[e]
        w = edge_weights[e] if edge_weights != NULL else 1.0

        if src < n_nodes:
            strengths[src] += w
        if dst < n_nodes and dst != src:
            strengths[dst] += w

        total_weight += w

    if total_weight < EPSILON:
        free(strengths)
        return 0.0

    cdef double m2 = total_weight  # 2m for undirected

    # Compute modularity
    cdef double Q = 0.0
    cdef double expected

    for e in range(n_edges):
        src = edge_src[e]
        dst = edge_dst[e]
        w = edge_weights[e] if edge_weights != NULL else 1.0

        if src < n_nodes and dst < n_nodes:
            # Only count if same community
            if partition[src] == partition[dst] and partition[src] >= 0:
                expected = (strengths[src] * strengths[dst]) / m2
                Q += w - expected

    free(strengths)

    return Q / m2


# ============================================================================
# Assortativity (categorical)
# ============================================================================

cdef double compute_assortativity_categorical(
    uint32_t* edge_src,
    uint32_t* edge_dst,
    double* edge_weights,
    uint64_t n_edges,
    int32_t* categories,
    uint32_t n_nodes,
) noexcept nogil:
    """
    Compute Newman's assortativity coefficient for categorical attribute.

    r = (Tr(e) - ||e²||) / (1 - ||e²||)

    where e_{ab} is fraction of edges connecting categories a and b.
    """
    if n_edges == 0 or n_nodes == 0:
        return 0.0

    # Find number of unique categories
    cdef int32_t max_cat = -1
    cdef uint32_t i

    for i in range(n_nodes):
        if categories[i] > max_cat:
            max_cat = categories[i]

    if max_cat < 0:
        return 0.0

    cdef uint32_t n_cats = <uint32_t>(max_cat + 1)

    # Build mixing matrix e_{ab}
    cdef double* mixing = <double*>malloc(n_cats * n_cats * sizeof(double))
    cdef double* row_sums = <double*>malloc(n_cats * sizeof(double))

    if mixing == NULL or row_sums == NULL:
        if mixing != NULL:
            free(mixing)
        if row_sums != NULL:
            free(row_sums)
        return 0.0

    memset(mixing, 0, n_cats * n_cats * sizeof(double))
    memset(row_sums, 0, n_cats * sizeof(double))

    cdef double total_weight = 0.0
    cdef uint64_t e
    cdef uint32_t src, dst
    cdef int32_t cat_src, cat_dst
    cdef double w

    for e in range(n_edges):
        src = edge_src[e]
        dst = edge_dst[e]
        w = edge_weights[e] if edge_weights != NULL else 1.0

        if src < n_nodes and dst < n_nodes:
            cat_src = categories[src]
            cat_dst = categories[dst]

            if cat_src >= 0 and cat_dst >= 0 and cat_src < <int32_t>n_cats and cat_dst < <int32_t>n_cats:
                mixing[cat_src * n_cats + cat_dst] += w
                if src != dst:
                    mixing[cat_dst * n_cats + cat_src] += w
                total_weight += w

    if total_weight < EPSILON:
        free(mixing)
        free(row_sums)
        return 0.0

    # Normalize mixing matrix
    cdef uint32_t a, b
    for a in range(n_cats):
        for b in range(n_cats):
            mixing[a * n_cats + b] /= total_weight

    # Compute row sums
    for a in range(n_cats):
        for b in range(n_cats):
            row_sums[a] += mixing[a * n_cats + b]

    # Compute Tr(e) and ||e²||
    cdef double trace = 0.0
    cdef double sum_sq = 0.0

    for a in range(n_cats):
        trace += mixing[a * n_cats + a]
        sum_sq += row_sums[a] * row_sums[a]

    free(mixing)
    free(row_sums)

    if fabs(1.0 - sum_sq) < EPSILON:
        return 0.0

    return (trace - sum_sq) / (1.0 - sum_sq)


# ============================================================================
# Neighbor taxonomy entropy
# ============================================================================

cdef double compute_neighbor_tax_entropy(
    int32_t* neighbor_taxa,
    double* neighbor_weights,
    uint32_t n_neighbors,
    uint32_t* n_unique_taxa,  # Output: number of unique taxa
) noexcept nogil:
    """
    Compute entropy of taxonomic labels among weighted neighbors.

    H_neighbor(j) = -Σ_a P_j(a) log P_j(a)

    where P_j(a) = (1/k_j) Σ_{k: T_k = a} A_{jk}
    """
    if n_neighbors == 0:
        if n_unique_taxa != NULL:
            n_unique_taxa[0] = 0
        return 0.0

    # Count unique taxa and compute weights per taxon
    # First pass: find max taxon ID
    cdef int32_t max_tax = -1
    cdef uint32_t i

    for i in range(n_neighbors):
        if neighbor_taxa[i] > max_tax:
            max_tax = neighbor_taxa[i]

    if max_tax < 0:
        if n_unique_taxa != NULL:
            n_unique_taxa[0] = 0
        return 0.0

    cdef uint32_t n_taxa = <uint32_t>(max_tax + 1)

    # Accumulate weights per taxon
    cdef double* tax_weights = <double*>malloc(n_taxa * sizeof(double))
    if tax_weights == NULL:
        if n_unique_taxa != NULL:
            n_unique_taxa[0] = 0
        return 0.0

    memset(tax_weights, 0, n_taxa * sizeof(double))

    cdef double total_weight = 0.0
    cdef double w
    cdef int32_t tax

    for i in range(n_neighbors):
        tax = neighbor_taxa[i]
        w = neighbor_weights[i] if neighbor_weights != NULL else 1.0
        if tax >= 0 and tax < <int32_t>n_taxa:
            tax_weights[tax] += w
            total_weight += w

    if total_weight < EPSILON:
        free(tax_weights)
        if n_unique_taxa != NULL:
            n_unique_taxa[0] = 0
        return 0.0

    # Count unique taxa and compute entropy
    cdef uint32_t unique_count = 0
    cdef double entropy = 0.0
    cdef double p

    for i in range(n_taxa):
        if tax_weights[i] > EPSILON:
            unique_count += 1
            p = tax_weights[i] / total_weight
            entropy -= p * log2(p)

    free(tax_weights)

    if n_unique_taxa != NULL:
        n_unique_taxa[0] = unique_count

    return entropy


# ============================================================================
# Dataset-level QC metrics
# ============================================================================

cdef NetworkQCDataset compute_dataset_qc_metrics(
    double* weights,
    uint32_t* read_counts,
    uint32_t n_refs,
    int32_t* community_ids,
    int32_t* taxonomy_ids,
    uint32_t* edge_list_src,
    uint32_t* edge_list_dst,
    double* edge_weights,
    uint64_t n_edges,
    NetworkQCConfig* config,
) noexcept nogil:
    """
    Compute all dataset-level QC metrics.
    """
    cdef NetworkQCDataset result
    memset(&result, 0, sizeof(NetworkQCDataset))

    if n_refs == 0:
        return result

    # Convert read counts to doubles for analysis
    cdef double* read_counts_d = <double*>malloc(n_refs * sizeof(double))
    if read_counts_d == NULL:
        return result

    cdef uint32_t i
    cdef double total_w = 0.0

    for i in range(n_refs):
        read_counts_d[i] = <double>read_counts[i] if read_counts != NULL else 0.0

    # Diversity metrics
    if weights != NULL:
        result.n_eff_weights = compute_n_eff(weights, n_refs)
        result.gini_weights = compute_gini(weights, n_refs)
        result.entropy_weights = compute_entropy(weights, n_refs)

        # Find max and compute top-k shares
        total_w = 0.0
        result.max_weight = 0.0
        for i in range(n_refs):
            if weights[i] > result.max_weight:
                result.max_weight = weights[i]
            total_w += weights[i]

        # For top-k, we'd need to sort - simplified here
        # (In full implementation, sort and sum top 5/10)

    result.n_eff_reads = compute_n_eff(read_counts_d, n_refs)
    result.gini_reads = compute_gini(read_counts_d, n_refs)

    free(read_counts_d)

    # Network metrics (if edges provided)
    if n_edges > 0 and edge_list_src != NULL and edge_list_dst != NULL:
        if taxonomy_ids != NULL:
            result.modularity_taxonomy = compute_modularity(
                edge_list_src, edge_list_dst, edge_weights, n_edges,
                taxonomy_ids, n_refs
            )
            result.assortativity_taxonomy = compute_assortativity_categorical(
                edge_list_src, edge_list_dst, edge_weights, n_edges,
                taxonomy_ids, n_refs
            )

        if community_ids != NULL:
            result.modularity_community = compute_modularity(
                edge_list_src, edge_list_dst, edge_weights, n_edges,
                community_ids, n_refs
            )

    return result


# ============================================================================
# Per-reference QC metrics
# ============================================================================

cdef void compute_reference_qc_metrics(
    uint32_t ref_idx,
    NetworkQCReference* result,
    int32_t* community_ids,
    int32_t* taxonomy_ids,
    uint32_t* neighbors,
    double* neighbor_weights,
    uint32_t n_neighbors,
    double* community_strength_sums,
    double* community_strength_sq,
    uint32_t* community_counts,
    uint32_t n_communities,
    int32_t ref_community,
    NetworkQCConfig* config,
) noexcept nogil:
    """
    Compute QC metrics for a single reference.

    Simplified version that only computes neighbor_tax_entropy and edge statistics.
    The tax_ambiguity_flag is computed separately by compute_tax_ambiguity_flags()
    which takes structural_role into account (HUBs are not flagged).
    """
    # All cdef declarations at top (required for nogil)
    cdef double total_strength = 0.0
    cdef uint32_t i
    cdef double w
    cdef double* comm_strengths = NULL
    cdef int32_t neigh_comm
    cdef double within_comm_strength
    cdef double mean_strength
    cdef double var
    cdef int32_t* neigh_taxa = NULL
    cdef uint32_t n_unique_taxa = 0

    memset(result, 0, sizeof(NetworkQCReference))

    if n_neighbors == 0:
        return

    # Compute total strength to neighbors
    for i in range(n_neighbors):
        w = neighbor_weights[i] if neighbor_weights != NULL else 1.0
        total_strength += w

    # Cross-community edge count and within-community z-score
    if community_ids != NULL and n_communities > 0:
        comm_strengths = <double*>malloc(n_communities * sizeof(double))
        if comm_strengths != NULL:
            memset(comm_strengths, 0, n_communities * sizeof(double))

            for i in range(n_neighbors):
                neigh_comm = community_ids[neighbors[i]]
                w = neighbor_weights[i] if neighbor_weights != NULL else 1.0
                if neigh_comm >= 0 and neigh_comm < <int32_t>n_communities:
                    comm_strengths[neigh_comm] += w
                    if neigh_comm != ref_community:
                        result.cross_community_edges += 1

            # Within-community z-score
            if ref_community >= 0 and ref_community < <int32_t>n_communities:
                within_comm_strength = comm_strengths[ref_community]
                if community_counts != NULL and community_counts[ref_community] > 1:
                    mean_strength = community_strength_sums[ref_community] / community_counts[ref_community]
                    var = (community_strength_sq[ref_community] / community_counts[ref_community]) - (mean_strength * mean_strength)
                    if var > EPSILON:
                        result.within_community_zscore = <float>((within_comm_strength - mean_strength) / sqrt(var))

            free(comm_strengths)

    # Neighbor taxonomic entropy (key metric for ambiguity detection)
    if taxonomy_ids != NULL:
        neigh_taxa = <int32_t*>malloc(n_neighbors * sizeof(int32_t))
        if neigh_taxa != NULL:
            for i in range(n_neighbors):
                neigh_taxa[i] = taxonomy_ids[neighbors[i]]

            n_unique_taxa = 0
            result.neighbor_tax_entropy = <float>compute_neighbor_tax_entropy(
                neigh_taxa, neighbor_weights, n_neighbors, &n_unique_taxa
            )

            free(neigh_taxa)

    # Edge weight statistics
    result.mean_edge_weight = <float>(total_strength / n_neighbors) if n_neighbors > 0 else 0.0
    result.max_edge_weight = 0.0
    for i in range(n_neighbors):
        w = neighbor_weights[i] if neighbor_weights != NULL else 1.0
        if w > result.max_edge_weight:
            result.max_edge_weight = <float>w


# ============================================================================
# Taxonomic ambiguity flagging
# ============================================================================

cdef void compute_tax_ambiguity_flags(
    NetworkQCReference* ref_metrics,
    uint32_t n_refs,
    char* structural_roles,
    NetworkQCConfig* config,
) noexcept nogil:
    """
    Compute taxonomic ambiguity flags based on neighbor entropy.

    This replaces the old chimera detection which was misnamed.
    HUB references (structural_role=2) are skipped since they are
    expected to have high entropy by nature.

    Flag levels based on neighbor_tax_entropy (in bits):
    - 0 = clean: entropy < biased_threshold (default: 1.0)
    - 1 = biased: entropy >= biased_threshold
    - 2 = mixed: entropy >= mixed_threshold (default: 2.0)
    - 3 = highly_mixed: entropy >= highly_mixed_threshold (default: 3.0)
    """
    if ref_metrics == NULL or n_refs == 0 or config == NULL:
        return

    cdef uint32_t i
    cdef double entropy
    cdef double biased_thresh = config.entropy_biased_threshold
    cdef double mixed_thresh = config.entropy_mixed_threshold
    cdef double highly_mixed_thresh = config.entropy_highly_mixed_threshold
    cdef char role

    # Use defaults if thresholds not set
    if biased_thresh <= 0:
        biased_thresh = 1.0
    if mixed_thresh <= 0:
        mixed_thresh = 2.0
    if highly_mixed_thresh <= 0:
        highly_mixed_thresh = 3.0

    for i in range(n_refs):
        # Skip HUB references (role=2) - they're expected to have high entropy
        if structural_roles != NULL:
            role = structural_roles[i]
            if role == 2:  # HUB
                ref_metrics[i].tax_ambiguity_flag = 0  # Clean (expected behavior)
                continue

        # Get raw neighbor entropy (in bits)
        entropy = ref_metrics[i].neighbor_tax_entropy

        # Set flag based on entropy thresholds
        if entropy >= highly_mixed_thresh:
            ref_metrics[i].tax_ambiguity_flag = 3  # Highly mixed
        elif entropy >= mixed_thresh:
            ref_metrics[i].tax_ambiguity_flag = 2  # Mixed
        elif entropy >= biased_thresh:
            ref_metrics[i].tax_ambiguity_flag = 1  # Biased
        else:
            ref_metrics[i].tax_ambiguity_flag = 0  # Clean


# ============================================================================
# Python-accessible functions
# ============================================================================

def compute_gini_py(values):
    """Python wrapper for Gini coefficient computation."""
    import numpy as np
    cdef double[::1] arr = np.ascontiguousarray(values, dtype=np.float64)
    return compute_gini(&arr[0], len(arr))


def compute_entropy_py(values):
    """Python wrapper for entropy computation."""
    import numpy as np
    cdef double[::1] arr = np.ascontiguousarray(values, dtype=np.float64)
    return compute_entropy(&arr[0], len(arr))


def compute_n_eff_py(values):
    """Python wrapper for effective number computation."""
    import numpy as np
    cdef double[::1] arr = np.ascontiguousarray(values, dtype=np.float64)
    return compute_n_eff(&arr[0], len(arr))


def compute_modularity_py(edge_src, edge_dst, edge_weights, partition, n_nodes):
    """Python wrapper for modularity computation."""
    import numpy as np
    cdef uint32_t[::1] src = np.ascontiguousarray(edge_src, dtype=np.uint32)
    cdef uint32_t[::1] dst = np.ascontiguousarray(edge_dst, dtype=np.uint32)
    cdef double[::1] weights = np.ascontiguousarray(edge_weights, dtype=np.float64)
    cdef int32_t[::1] part = np.ascontiguousarray(partition, dtype=np.int32)

    return compute_modularity(
        &src[0], &dst[0], &weights[0], len(src),
        &part[0], n_nodes
    )


def compute_assortativity_py(edge_src, edge_dst, edge_weights, categories, n_nodes):
    """Python wrapper for assortativity computation."""
    import numpy as np
    cdef uint32_t[::1] src = np.ascontiguousarray(edge_src, dtype=np.uint32)
    cdef uint32_t[::1] dst = np.ascontiguousarray(edge_dst, dtype=np.uint32)
    cdef double[::1] weights = np.ascontiguousarray(edge_weights, dtype=np.float64)
    cdef int32_t[::1] cats = np.ascontiguousarray(categories, dtype=np.int32)

    return compute_assortativity_categorical(
        &src[0], &dst[0], &weights[0], len(src),
        &cats[0], n_nodes
    )
