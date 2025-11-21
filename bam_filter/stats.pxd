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
    void* trusted_reads_hash_int
) nogil
