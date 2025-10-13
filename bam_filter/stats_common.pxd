# cython: language_level=3
from libc.stdint cimport int32_t, int64_t

# Shared C-level definitions used across stats modules
# Keep only minimal, widely-used structs here to avoid circular imports.

cdef struct ref_stats_t:
    # Basic counts
    int64_t n_reads
    int64_t n_alns

    # Read statistics
    double read_length_mean
    double read_length_std
    int read_length_median
    int min_read_length
    int max_read_length
    int read_length_mode

    # ANI and quality metrics
    double ani_mean
    double ani_std
    double ani_median
    double min_ani
    double max_ani
    double aligned_length_mean
    double aln_score_mean
    double aln_score_std
    double mapq_mean
    double mapq_std
    double edit_dist_mean
    double edit_dist_std
    double gc_content_mean
    double gc_content_std
    double gc_content_total

    # Coverage statistics
    int64_t bases_covered
    int64_t total_coverage
    double mean_coverage
    double mean_coverage_covered
    double breadth
    double exp_breadth
    double breadth_exp_ratio
    double cov_evenness
    double c_v
    double d_i

    # TAD (Truncated Average Depth) stats
    double mean_coverage_trunc
    int64_t mean_coverage_trunc_len
    int64_t n_reads_tad
    int64_t tax_abund_aln
    int64_t tax_abund_read
    int64_t tax_abund_tad

    # Coverage distribution stats
    double entropy
    double norm_entropy
    double gini
    double norm_gini
    int64_t n_bins
    double site_density

    # Interval merging results
    int64_t max_covered_bases
    double mean_covered_bases

    # Reference lengths
    int64_t ref_length
    int64_t bam_ref_length
