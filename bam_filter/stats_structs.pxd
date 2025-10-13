# Shared C struct declarations for stats modules
from libc.stdint cimport int64_t

cdef struct ref_stats_t:
    int64_t n_reads
    int64_t n_alns
    double read_length_mean
    double read_length_std
    int read_length_median
    int min_read_length
    int max_read_length
    int read_length_mode
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
    double mean_coverage_trunc
    int64_t mean_coverage_trunc_len
    int64_t n_reads_tad
    int64_t tax_abund_aln
    int64_t tax_abund_read
    int64_t tax_abund_tad
    double entropy
    double norm_entropy
    double gini
    double norm_gini
    int64_t n_bins
    double site_density
    int64_t max_covered_bases
    double mean_covered_bases
    int64_t ref_length
    int64_t bam_ref_length


cdef struct filter_conditions_t:
    int min_read_count
    double min_avg_read_ani
    double min_expected_breadth_ratio
    double min_breadth
    double min_coverage_evenness
    double max_coeff_var
    double min_coverage_mean
    double min_norm_entropy
    double max_norm_gini
    bint enable_min_avg_read_ani
    bint enable_min_expected_breadth_ratio
    bint enable_min_breadth
    bint enable_min_coverage_evenness
    bint enable_max_coeff_var
    bint enable_min_coverage_mean
    bint enable_min_norm_entropy
    bint enable_max_norm_gini
