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
Taxonomy-informed filtering logic.

This module implements multi-evidence filtering that combines:
1. Graph topology (clustering coefficient, betweenness centrality)
2. Community structure (Community communities)
3. Taxonomic congruence (taxonomy flags)
4. Statistical outlier detection

FILTERING STRATEGY:
- Level 1 (Strict): Automatic removal of high-confidence contamination
- Level 2 (Informed): Taxonomy-aware outlier detection
- Level 3 (Validation): Second-chance evaluation for graph-flagged references
"""

from libc.stdio cimport printf
from libc.stdint cimport uint32_t, int32_t, uint8_t
from libc.math cimport fabs

from bam_filter.processor_graph cimport ReferencePattern

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil


# ==============================================================================
# Configuration Structure
# ==============================================================================

cdef struct TaxonomyFilterConfig:
    # Enable taxonomy-informed filtering
    bint enabled

    # Level 1: Strict filtering thresholds
    bint enable_strict_filtering              # Enable automatic removal
    uint32_t strict_min_connections           # Min connections for strict filtering (default: 5)
    float strict_cross_domain_fraction        # Cross-domain fraction threshold (default: 0.10)

    # Level 2: Taxonomy-aware outlier detection
    bint enable_weighted_outlier_detection    # Weight outlier scores by taxonomy
    float taxonomy_anomaly_weight             # Weight multiplier for taxonomy flags (default: 2.0)

    # Level 3: Second-chance validation
    bint enable_second_chance                 # Give flagged refs with normal taxonomy another chance
    float second_chance_cc_threshold          # Min CC for second chance (default: 0.3)


# ==============================================================================
# Level 1: Strict Filtering (High-Confidence Contamination)
# ==============================================================================

cdef uint32_t apply_strict_taxonomy_filter(
    ReferencePattern* pattern_data,
    uint32_t n_refs,
    char* keep_flag,
    TaxonomyFilterConfig* config,
    uint32_t* connection_counts,
    bint verbose
) noexcept nogil:
    """
    Apply strict taxonomy-based filtering for high-confidence contamination.

    CRITERIA for automatic removal:
    1. Cross-domain contamination (taxonomy_flag == 2)
    2. High connectivity (>= strict_min_connections)
    3. Not already removed by Community

    This catches obvious contamination cases where a reference:
    - Connects organisms from different domains (Bacteria <-> Eukaryota)
    - Has many graph connections (not just one spurious alignment)

    Parameters
    ----------
    pattern_data : ReferencePattern*
        Array of reference patterns with taxonomy info
    n_refs : uint32_t
        Number of references
    keep_flag : char*
        Keep flags (1=keep, 0=remove), modified in-place
    config : TaxonomyFilterConfig*
        Configuration parameters
    connection_counts : uint32_t*
        Per-reference connection counts from graph analysis
    verbose : bint
        Enable verbose logging

    Returns
    -------
    uint32_t
        Number of references removed by strict filtering
    """
    if not config.enabled or not config.enable_strict_filtering:
        return 0

    cdef uint32_t ref_idx
    cdef uint32_t removed_count = 0
    cdef uint32_t cross_domain_count = 0
    cdef uint32_t kingdom_mismatch_count = 0
    cdef char tax_flag
    cdef uint32_t connections

    if verbose:
        bf_nogil_logf_notime(
            b"TAXONOMY_FILTER",
            "\n=== STRICT TAXONOMY FILTERING ===\n"
        )
        bf_nogil_logf_notime(
            b"TAXONOMY_FILTER",
            "  Min connections for strict filtering: %u\n",
            config.strict_min_connections
        )

    for ref_idx in range(n_refs):
        # Skip if already removed by Community
        if keep_flag[ref_idx] == 0:
            continue

        tax_flag = pattern_data[ref_idx].taxonomy_flag
        connections = 0
        if connection_counts != NULL:
            connections = connection_counts[ref_idx]

        # Check for cross-domain contamination with high connectivity
        if tax_flag == 2:  # cross_domain
            cross_domain_count += 1
            if connections >= config.strict_min_connections:
                keep_flag[ref_idx] = 0
                removed_count += 1

        # Optionally also remove kingdom mismatches with very high connectivity
        elif tax_flag == 3:  # kingdom_mismatch
            kingdom_mismatch_count += 1
            # More conservative: require more connections for kingdom mismatch
            if connections >= config.strict_min_connections * 2:
                keep_flag[ref_idx] = 0
                removed_count += 1

    if verbose:
        bf_nogil_logf_notime(
            b"TAXONOMY_FILTER",
            "  Found %u cross-domain flagged references\n",
            cross_domain_count
        )
        bf_nogil_logf_notime(
            b"TAXONOMY_FILTER",
            "  Found %u kingdom-mismatch flagged references\n",
            kingdom_mismatch_count
        )
        bf_nogil_logf_notime(
            b"TAXONOMY_FILTER",
            "  REMOVED %u references by strict taxonomy filtering\n",
            removed_count
        )

    return removed_count


# ==============================================================================
# Level 2: Taxonomy-Aware Anomaly Score Weighting
# ==============================================================================

cdef uint32_t weight_anomaly_scores_by_taxonomy(
    ReferencePattern* pattern_data,
    uint32_t n_refs,
    float* anomaly_scores,
    TaxonomyFilterConfig* config,
    bint verbose
) noexcept nogil:
    """
    Weight anomaly scores by taxonomy flags to increase sensitivity.

    RATIONALE:
    - References with taxonomy anomalies should be more suspicious
    - Increase their anomaly scores so outlier detection is more likely to flag them
    - This creates a synergy between graph and taxonomy evidence

    WEIGHTING SCHEME:
    - taxonomy_flag == 0 (normal): no change
    - taxonomy_flag == 1 (genus mismatch): multiply by taxonomy_anomaly_weight
    - taxonomy_flag == 2 (cross-domain): multiply by taxonomy_anomaly_weight * 2
    - taxonomy_flag == 3 (kingdom mismatch): multiply by taxonomy_anomaly_weight * 1.5

    Parameters
    ----------
    pattern_data : ReferencePattern*
        Array of reference patterns with taxonomy info
    n_refs : uint32_t
        Number of references
    anomaly_scores : float*
        Anomaly scores from outlier detection, modified in-place
    config : TaxonomyFilterConfig*
        Configuration parameters
    verbose : bint
        Enable verbose logging
    """
    if not config.enabled or not config.enable_weighted_outlier_detection:
        return 0

    if anomaly_scores == NULL:
        return 0

    cdef uint32_t ref_idx
    cdef char tax_flag
    cdef float weight_multiplier
    cdef float original_score
    cdef uint32_t weighted_count = 0

    if verbose:
        bf_nogil_logf_notime(
            b"TAXONOMY_FILTER",
            "\n=== TAXONOMY-WEIGHTED ANOMALY SCORING ===\n"
        )
        bf_nogil_logf_notime(
            b"TAXONOMY_FILTER",
            "  Anomaly weight multiplier: %.2f\n",
            config.taxonomy_anomaly_weight
        )

    for ref_idx in range(n_refs):
        tax_flag = pattern_data[ref_idx].taxonomy_flag

        if tax_flag == 0:
            continue  # Normal taxonomy, no weighting

        original_score = anomaly_scores[ref_idx]

        # Determine weight multiplier based on taxonomy flag severity
        if tax_flag == 2:  # cross_domain (most severe)
            weight_multiplier = config.taxonomy_anomaly_weight * 2.0
        elif tax_flag == 3:  # kingdom_mismatch
            weight_multiplier = config.taxonomy_anomaly_weight * 1.5
        elif tax_flag == 1:  # genus_mismatch (least severe)
            weight_multiplier = config.taxonomy_anomaly_weight
        else:
            weight_multiplier = 1.0

        # Apply weight (ensure score stays in valid range [0, 1])
        anomaly_scores[ref_idx] = original_score * weight_multiplier
        if anomaly_scores[ref_idx] > 1.0:
            anomaly_scores[ref_idx] = 1.0

        weighted_count += 1

    if verbose:
        bf_nogil_logf_notime(
            b"TAXONOMY_FILTER",
            "  Weighted %u references with taxonomy anomalies\n",
            weighted_count
        )

    return weighted_count


# ==============================================================================
# Level 3: Second-Chance Validation
# ==============================================================================

cdef uint32_t apply_second_chance_validation(
    ReferencePattern* pattern_data,
    uint32_t n_refs,
    char* keep_flag,
    TaxonomyFilterConfig* config,
    bint verbose
) noexcept nogil:
    """
    Give references flagged by graph analysis a second chance if they have normal taxonomy.

    RATIONALE:
    - Graph-based filtering (Community clustering) can have false positives
    - References with low CC but NORMAL taxonomy might be legitimate
    - Especially important for:
      - Conserved genes that naturally have hub-like patterns
      - High-copy sequences (rRNA, transposons)
      - References with genuinely diverse read mappings

    CRITERIA for second chance (restore to keep_flag=1):
    1. Currently flagged for removal (keep_flag == 0)
    2. Normal taxonomy (taxonomy_flag == 0)
    3. Individual CC above threshold (community_individual_cc >= second_chance_cc_threshold)
    4. Not removed by strict filtering (taxonomy has to be clean)

    Parameters
    ----------
    pattern_data : ReferencePattern*
        Array of reference patterns with taxonomy and Community info
    n_refs : uint32_t
        Number of references
    keep_flag : char*
        Keep flags (1=keep, 0=remove), modified in-place
    config : TaxonomyFilterConfig*
        Configuration parameters
    verbose : bint
        Enable verbose logging

    Returns
    -------
    uint32_t
        Number of references restored by second-chance validation
    """
    if not config.enabled or not config.enable_second_chance:
        return 0

    cdef uint32_t ref_idx
    cdef uint32_t restored_count = 0
    cdef char tax_flag
    cdef float individual_cc

    if verbose:
        bf_nogil_logf_notime(
            b"TAXONOMY_FILTER",
            "\n=== SECOND-CHANCE VALIDATION ===\n"
        )
        bf_nogil_logf_notime(
            b"TAXONOMY_FILTER",
            "  Min CC for second chance: %.3f\n",
            config.second_chance_cc_threshold
        )

    for ref_idx in range(n_refs):
        # Only consider references flagged for removal
        if keep_flag[ref_idx] == 1:
            continue

        tax_flag = pattern_data[ref_idx].taxonomy_flag
        individual_cc = pattern_data[ref_idx].community_individual_cc

        # Restore if: normal taxonomy AND decent clustering coefficient
        if tax_flag == 0 and individual_cc >= config.second_chance_cc_threshold:
            keep_flag[ref_idx] = 1
            restored_count += 1

    if verbose:
        bf_nogil_logf_notime(
            b"TAXONOMY_FILTER",
            "  RESTORED %u references with normal taxonomy and CC >= %.3f\n",
            restored_count,
            config.second_chance_cc_threshold
        )

    return restored_count


# ==============================================================================
# Master Function: Apply All Taxonomy-Informed Filtering
# ==============================================================================

cdef TaxonomyFilterStats apply_taxonomy_informed_filtering(
    ReferencePattern* pattern_data,
    uint32_t n_refs,
    char* keep_flag,
    float* anomaly_scores,
    uint32_t* connection_counts,
    TaxonomyFilterConfig* config,
    bint verbose
) noexcept nogil:
    """
    Apply all levels of taxonomy-informed filtering.

    ORDER OF OPERATIONS:
    1. Level 2: Weight anomaly scores by taxonomy (before Community makes decisions)
    2. (Community clustering runs with weighted scores)
    3. Level 1: Strict filtering for high-confidence contamination
    4. Level 3: Second-chance validation for false positives

    Parameters
    ----------
    pattern_data : ReferencePattern*
        Array of reference patterns with taxonomy and graph info
    n_refs : uint32_t
        Number of references
    keep_flag : char*
        Keep flags (1=keep, 0=remove), modified in-place by Community and this function
    anomaly_scores : float*
        Anomaly scores from outlier detection (will be weighted if enabled)
    connection_counts : uint32_t*
        Per-reference connection counts from graph analysis
    config : TaxonomyFilterConfig*
        Configuration parameters
    verbose : bint
        Enable verbose logging
    """
    cdef TaxonomyFilterStats stats
    stats.strict_removed = 0
    stats.weighted_count = 0
    stats.second_chance_restored = 0

    if not config.enabled:
        if verbose:
            bf_nogil_logf_notime(
                b"TAXONOMY_FILTER",
                "Taxonomy-informed filtering disabled\n"
            )
        return stats

    cdef uint32_t strict_removed = 0
    cdef uint32_t second_chance_restored = 0

    if verbose:
        bf_nogil_logf_notime(
            b"TAXONOMY_FILTER",
            "\n================================================\n"
        )
        bf_nogil_logf_notime(
            b"TAXONOMY_FILTER",
            "TAXONOMY-INFORMED FILTERING\n"
        )
        bf_nogil_logf_notime(
            b"TAXONOMY_FILTER",
            "================================================\n"
        )

    # Note: Level 2 (weight_anomaly_scores_by_taxonomy) should be called BEFORE
    # Community clustering runs, so it's exposed separately. Here we only do
    # post-Community filtering (Levels 1 and 3).

    # Level 1: Strict filtering (remove high-confidence contamination)
    strict_removed = apply_strict_taxonomy_filter(
        pattern_data, n_refs, keep_flag, config, connection_counts, verbose
    )

    # Level 3: Second-chance validation (restore false positives)
    second_chance_restored = apply_second_chance_validation(
        pattern_data, n_refs, keep_flag, config, verbose
    )

    # Populate stats structure
    stats.strict_removed = strict_removed
    stats.second_chance_restored = second_chance_restored
    # weighted_count is set by weight_anomaly_scores_by_taxonomy (called before this function)

    if verbose:
        bf_nogil_logf_notime(
            b"TAXONOMY_FILTER",
            "\n=== TAXONOMY FILTERING SUMMARY ===\n"
        )
        bf_nogil_logf_notime(
            b"TAXONOMY_FILTER",
            "  Strict filtering removed: %u\n",
            strict_removed
        )
        bf_nogil_logf_notime(
            b"TAXONOMY_FILTER",
            "  Second-chance restored: %u\n",
            second_chance_restored
        )
        bf_nogil_logf_notime(
            b"TAXONOMY_FILTER",
            "  Net impact: %d references\n",
            <int32_t>second_chance_restored - <int32_t>strict_removed
        )
        bf_nogil_logf_notime(
            b"TAXONOMY_FILTER",
            "================================================\n\n"
        )

    return stats


# ==============================================================================
# Python-Level Interface
# ==============================================================================

def py_apply_taxonomy_filtering(
    pattern_data,  # Would need proper Python wrapper type
    keep_flag_array,
    config_dict,
    verbose=False
):
    """
    Python wrapper for taxonomy-informed filtering (for testing).

    This is a placeholder - actual integration happens at C level.
    """
    raise NotImplementedError("Use C-level integration in processor.pyx")
