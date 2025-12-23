# cython: language_level=3
from libc.stdint cimport int32_t, int64_t, uint32_t
from bam_filter.processor_types cimport samFile, sam_hdr_t, hts_idx_t

# RLE coverage data structures
cdef struct RLEInterval:
    int64_t start
    int64_t end
    int32_t count

cdef struct RLECoverage:
    RLEInterval* intervals
    int64_t n_intervals
    int64_t capacity
    int64_t ref_length

cdef extern from "seqid_khash.h":
    ctypedef struct kh_seqid_map_t

cdef extern from "taxonomy_khash.h":
    ctypedef uint32_t khint_t
    ctypedef struct kh_str_t:
        pass


# RefStats structure (declared here so callers can pass pointers to C functions)
cdef struct RefStats:
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

    # ANI and quality metrics (raw, from edit distance)
    double ani_mean
    double ani_std
    double ani_median
    double min_ani
    double max_ani

    # Corrected ANI (damage-corrected, set by PMD stage)
    double ani_corrected_mean
    double ani_corrected_std

    # Per-reference damage counts for computing corrected ANI (stored in Pass 1)
    int64_t total_aligned_length    # Sum of aligned_length across all alignments
    int64_t total_match_count       # Sum of match_count across all alignments
    int64_t total_ct_5p_count       # Sum of C→T mismatches in first 8bp (5' end)
    int64_t total_ga_3p_count       # Sum of G→A mismatches in last 8bp (3' end)
    int64_t total_c_at_5p_count     # Sum of C bases in reference at first 8bp (damage opportunities)
    int64_t total_g_at_3p_count     # Sum of G bases in reference at last 8bp (damage opportunities)
    double sum_damage_llr           # Sum of position-specific damage log-likelihood ratios
    int64_t damage_llr_count        # Number of alignments with damage_llr computed
    # Per-position damage counts for full damage model (20 positions from each end)
    double n_5p[20]                 # C bases at each 5' position (opportunities)
    double k_5p[20]                 # C→T mismatches at each 5' position
    double n_3p[20]                 # G bases at each 3' position (opportunities)
    double k_3p[20]                 # G→A mismatches at each 3' position

    double aligned_length_mean
    double aln_score_mean
    double aln_score_std

    # ZP/ZS tag statistics (from reassign/EM output)
    double zp_mean     # Mean EM posterior probability
    double zp_std      # Std of EM posterior
    double zs_mean     # Mean log-likelihood alignment score
    double zs_std      # Std of log-likelihood score
    int64_t zp_count   # Number of alignments with ZP tag
    int64_t zs_count   # Number of alignments with ZS tag

    double mapq_mean
    double mapq_std
    double edit_dist_mean
    double edit_dist_std
    double read_gc_content_mean
    double read_gc_content_std
    double read_gc_content_total
    double ref_gc_content_mean
    double ref_gc_content_std
    double ref_gc_content_total
    double dust_mean
    double dust_std

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
    # Spatial entropy: measures uniformity of covered position distribution across reference
    # (NOT depth distribution - measures how evenly spread coverage is spatially)
    double spatial_entropy
    double norm_spatial_entropy
    double gini
    double norm_gini
    int64_t n_bins
    double site_density
    # Histogram summary cached from calculate_rle_coverage_stats
    int32_t hist_min
    int32_t hist_max
    int64_t hist_nonzero
    double hist_mean
    double hist_sd

    # Per-reference processing time (seconds)
    double ref_seconds
    # Breakdown of per-reference timing (seconds)
    double cov_events_sec
    double cov_merge_sec
    double cov_tad_sec
    double cov_total_sec
    double abundance_seconds
    # NOTE: histogram timing/debug fields removed (no per-ref histogram timing collected)

    # Interval merging results
    int64_t max_covered_bases
    double mean_covered_bases
    int64_t n_intervals  # Number of coverage intervals
    double sum_interval_length_sq  # Sum of interval_length^2 for WCB (computed inline)

    # Contamination detection metrics (computed post-hoc from existing stats)
    double weighted_contiguity_breadth  # WCB = sum(interval_len^2) / bases_covered^2 (Herfindahl index)
    double complexity_penalized_coverage  # CPC = breadth * (1 - dust_mean)
    double overlap_redundancy_index  # ORI = total_aligned_bases / bases_covered
    double mega_genome_sparsity_index  # MGSI = log10(expected_breadth / observed_breadth)
    double coverage_compressibility_ratio  # CCR = n_intervals / (bases_covered / read_len_mean)
    double feature_space_clustering_score  # FSCS = variance of GC/complexity across aligned regions

    # Authenticity metrics (computed post-hoc from spatial distribution)
    double authenticity_score  # norm_spatial_entropy - norm_gini (higher = more authentic)
    double authenticity_pvalue  # P(score <= x | real distribution), lower = likely contamination

    # Reference lengths
    int64_t ref_length
    int64_t bam_ref_length


cdef struct FilterConditions:
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

cdef void initialize_reference_stats(RefStats* stats, int64_t ref_length, int64_t bam_ref_length) noexcept nogil

# Import PMD types for damage correction
from bam_filter.processor_pmd cimport PMDStatsAccumulator, PMDCurve

cdef int calculate_reference_stats(
    samFile* htsfile,
    sam_hdr_t* header,
    hts_idx_t* idx,
    int64_t tid,
    int64_t num_alns,
    RefStats* stats,
    kh_seqid_map_t* unique_reads_map,
    double min_read_ani_c,
    int min_read_length_c,
    int max_read_length_c,
    int64_t scale,
    int trim_ends,
    int trim_min,
    int trim_max,
    bint verbose,
    void* trusted_reads_hash_int,
    PMDStatsAccumulator* pmd_acc,
    bint collect_damage_stats,
    PMDCurve* pmd_curve,
    float pmd_epsilon
) nogil
