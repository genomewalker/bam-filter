# cython: language_level=3
"""
C-level declarations for LCA stats processing.

This module calculates per-taxid statistics by aggregating coverage and quality
metrics across all references that contribute to each taxid after LCA assignment.
"""

from libc.stdint cimport int32_t, int64_t
from bam_filter.stats cimport RefStats, RLECoverage
from bam_filter.taxonomy_db cimport TaxonomyDB, TaxonomyDatabase

# Forward declarations
cdef struct TaxidInventory
cdef struct TaxidRefRLE
cdef struct TaxidBatchRLE
cdef struct PerRefStats
cdef struct TaxidAggregatedStats

# Data structures
# Note: kh_str_map_t is defined in seqid_khash.h and declared in processor_lca_stats.pyx
ctypedef void* kh_str_map_ptr  # Opaque pointer to kh_str_map_t
cdef struct TaxidInventory:
    int32_t taxid
    int32_t rank_id
    int64_t n_reads
    int64_t total_ref_length
    int32_t* ref_indices
    int32_t n_refs
    int32_t refs_capacity
    kh_str_map_ptr read_names  # Opaque pointer to kh_str_map_t hash

cdef struct TaxidRefRLE:
    int32_t ref_index
    int64_t ref_length
    RLECoverage* rle
    int64_t n_alns          # Total alignments to this reference
    int64_t read_length_sum # Sum of read lengths (for calculating mean)
    kh_str_map_ptr unique_reads  # Unique read names (for n_reads)

cdef struct TaxidBatchRLE:
    int32_t taxid
    TaxidRefRLE* ref_rles  # Contiguous array (not array of pointers!)
    int32_t n_refs
    int32_t capacity

cdef struct PerRefStats:
    int32_t ref_index
    char* ref_name
    int64_t n_reads
    RefStats stats

cdef struct TaxidAggregatedStats:
    int32_t taxid
    char* taxid_name
    int32_t rank_id
    char* rank_name
    int32_t n_refs
    int64_t total_reads
    int64_t total_alns
    int64_t total_bases

    # Read-weighted means
    double mean_read_length
    double mean_read_length_std
    double mean_read_length_median
    double mean_gc_content
    double mean_gc_content_std
    double mean_aligned_length
    double mean_aln_score
    double mean_mapq
    double mean_edit_dist
    double mean_ani
    double mean_ani_std
    double mean_ani_median
    int64_t bases_covered
    double mean_coverage
    double mean_tad_cov
    double mean_coverage_covered
    double mean_breadth
    double mean_exp_breadth
    double mean_breadth_exp_ratio
    double mean_entropy
    double mean_norm_entropy
    double mean_gini
    double mean_norm_gini
    double mean_cv
    double mean_di
    double mean_cov_evenness
    int64_t tax_abund_read
    int64_t tax_abund_aln
    int64_t tax_abund_tad
    int64_t n_reads_tad

    # Consistency metrics (CV across refs)
    double cv_coverage
    double cv_tad
    double cv_breadth
    double cv_breadth_exp_ratio
    double cv_entropy
    double cv_norm_entropy
    double cv_gini
    double cv_norm_gini
    double cv_cov_evenness
    double cv_ani
    double cv_gc_content
    double cv_aln_score
    double cv_mapq
    double cv_edit_dist
    double cv_ref_length

    # Reference length stats
    int64_t total_ref_length
    int64_t mean_ref_length
    int64_t min_ref_length
    int64_t max_ref_length

    # Best reference
    int32_t best_ref_index
    char* best_ref_name
    double best_coverage
    double best_breadth
    double best_tad
    double best_ani

    # Lineage string
    char* tax_path

    # Array of per-ref stats (for detailed output)
    PerRefStats* per_ref_stats

# Core functions
cdef int process_lca_stats(
    const char* bam_path,
    const char* lca_per_read_path,
    const char* output_path,
    TaxonomyDatabase taxdb_py,
    const char* taxonomy_db_path,
    int num_threads,
    bint verbose
) except -1
