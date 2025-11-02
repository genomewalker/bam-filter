# cython: language_level=3

"""
Taxonomy-informed filtering header.
"""

from libc.stdint cimport uint32_t, int32_t, uint8_t
from bam_filter.processor_graph cimport ReferencePattern


# Configuration structure for taxonomy-informed filtering
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


# Level 1: Strict filtering
cdef uint32_t apply_strict_taxonomy_filter(
    ReferencePattern* pattern_data,
    uint32_t n_refs,
    char* keep_flag,
    TaxonomyFilterConfig* config,
    uint32_t* connection_counts,
    bint verbose
) noexcept nogil


# Level 2: Anomaly score weighting
cdef uint32_t weight_anomaly_scores_by_taxonomy(
    ReferencePattern* pattern_data,
    uint32_t n_refs,
    float* anomaly_scores,
    TaxonomyFilterConfig* config,
    bint verbose
) noexcept nogil


# Level 3: Second-chance validation
cdef uint32_t apply_second_chance_validation(
    ReferencePattern* pattern_data,
    uint32_t n_refs,
    char* keep_flag,
    TaxonomyFilterConfig* config,
    bint verbose
) noexcept nogil


# Statistics structure for taxonomy filtering results
cdef struct TaxonomyFilterStats:
    uint32_t strict_removed
    uint32_t weighted_count
    uint32_t second_chance_restored


# Master function
cdef TaxonomyFilterStats apply_taxonomy_informed_filtering(
    ReferencePattern* pattern_data,
    uint32_t n_refs,
    char* keep_flag,
    float* anomaly_scores,
    uint32_t* connection_counts,
    TaxonomyFilterConfig* config,
    bint verbose
) noexcept nogil
