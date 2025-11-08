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

# IO implementation extracted from stats.pyx: TSV/BAM stats output writers

from libc.stdio cimport fprintf, FILE, fopen, fclose
from libc.string cimport strlen
from bam_filter.stats cimport RefStats, FilterConditions
from bam_filter.stats_bam_writer cimport passes_filters
from bam_filter.processor cimport sam_hdr_t, sam_hdr_tid2name


cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_verbose(int level, const char* tag, const char* fmt, ...) nogil
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil

# Fast endswith helper for .gz usable nogil (self-contained in IO module)
cdef inline bint c_str_endswith_gz(const char* s) nogil:
    if s == NULL:
        return False
    cdef int L = <int>strlen(s)
    if L < 3:
        return False
    # Check last three chars are '.' 'g' 'z'
    return s[L-3] == 46 and s[L-2] == 103 and s[L-1] == 122

cdef const char* UNKNOWN_REF_NAME = "<unknown>"

cdef inline const char* get_ref_name(sam_hdr_t* header, int tid) noexcept nogil:
    cdef const char* name = sam_hdr_tid2name(header, tid)
    if name == NULL:
        return UNKNOWN_REF_NAME
    return name

cdef int write_stats_to_file(
    FILE* fp,
    RefStats* global_ref_stats,
    sam_hdr_t* header,
    int n_refs,
    FilterConditions* filters,
    bint apply_filters
) except -1 nogil:
    """Helper function to write statistics to a file with optional filtering."""
    cdef int i
    cdef RefStats* stats
    cdef int refs_written = 0
    cdef const char* ref_name

    cdef const char* OUTPUT_HEADER = (
        "reference\t"
        "n_reads\t"
        "n_alns\t"
        "read_length_mean\t"
        "read_length_std\t"
        "read_length_min\t"
        "read_length_max\t"
        "read_length_median\t"
        "read_length_mode\t"
        "gc_content_mean\t"
        "gc_content_std\t"
        "gc_content_total\t"
        "dust_mean\t"
        "dust_std\t"
        "read_aligned_length\t"
        "read_aln_score\t"
        "mapping_quality\t"
        "edit_distances\t"
        "read_ani_mean\t"
        "read_ani_std\t"
        "read_ani_median\t"
        "bases_covered\t"
        "max_covered_bases\t"
        "mean_covered_bases\t"
        "coverage_mean\t"
        "coverage_mean_trunc\t"
        "coverage_mean_trunc_len\t"
        "coverage_covered_mean\t"
        "reference_length\t"
        "bam_reference_length\t"
        "breadth\t"
        "exp_breadth\t"
        "breadth_exp_ratio\t"
        "n_bins\t"
        "site_density\t"
        "entropy\t"
        "norm_entropy\t"
        "gini\t"
        "norm_gini\t"
        "c_v\t"
        "d_i\t"
        "cov_evenness\t"
        "tax_abund_read\t"
        "tax_abund_aln\t"
        "tax_abund_tad\t"
        "n_reads_tad\n"
    )

    cdef const char* OUTPUT_FORMAT = (
        "%s\t"
        "%lld\t"
        "%lld\t"
        "%.4f\t"
        "%.4f\t"
        "%d\t"
        "%d\t"
        "%d\t"
        "%d\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%lld\t"
        "%lld\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%lld\t"
        "%.4f\t"
        "%lld\t"
        "%lld\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%lld\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%lld\t"
        "%lld\t"
        "%lld\t"
        "%lld\n"
    )

    fprintf(fp, OUTPUT_HEADER)

    cdef int filtered_total = 0
    cdef int filtered_passed = 0

    for i in range(n_refs):
        stats = &global_ref_stats[i]
        if stats.n_alns > 0:
            filtered_total += 1
            ref_name = get_ref_name(header, i)
            if not apply_filters:
                fprintf(
                    fp,
                    OUTPUT_FORMAT,
                    ref_name,
                    <long long>stats.n_reads,
                    <long long>stats.n_alns,
                    stats.read_length_mean,
                    stats.read_length_std,
                    stats.min_read_length,
                    stats.max_read_length,
                    stats.read_length_median,
                    stats.read_length_mode,
                    stats.gc_content_mean,
                    stats.gc_content_std,
                    stats.gc_content_total,
                    stats.dust_mean,
                    stats.dust_std,
                    stats.aligned_length_mean,
                    stats.aln_score_mean,
                    stats.mapq_mean,
                    stats.edit_dist_mean,
                    stats.ani_mean,
                    stats.ani_std,
                    stats.ani_median,
                    <long long>stats.bases_covered,
                    <long long>stats.max_covered_bases,
                    stats.mean_covered_bases,
                    stats.mean_coverage,
                    stats.mean_coverage_trunc,
                    <long long>stats.mean_coverage_trunc_len,
                    stats.mean_coverage_covered,
                    <long long>stats.ref_length,
                    <long long>stats.bam_ref_length,
                    stats.breadth,
                    stats.exp_breadth,
                    stats.breadth_exp_ratio,
                    <long long>stats.n_bins,
                    stats.site_density,
                    stats.entropy,
                    stats.norm_entropy,
                    stats.gini,
                    stats.norm_gini,
                    stats.c_v,
                    stats.d_i,
                    stats.cov_evenness,
                    <long long>stats.tax_abund_read,
                    <long long>stats.tax_abund_aln,
                    <long long>stats.tax_abund_tad,
                    <long long>stats.n_reads_tad
                )
                refs_written += 1
                filtered_passed += 1
            elif passes_filters(stats, filters):
                fprintf(
                    fp,
                    OUTPUT_FORMAT,
                    ref_name,
                    <long long>stats.n_reads,
                    <long long>stats.n_alns,
                    stats.read_length_mean,
                    stats.read_length_std,
                    stats.min_read_length,
                    stats.max_read_length,
                    stats.read_length_median,
                    stats.read_length_mode,
                    stats.gc_content_mean,
                    stats.gc_content_std,
                    stats.gc_content_total,
                    stats.dust_mean,
                    stats.dust_std,
                    stats.aligned_length_mean,
                    stats.aln_score_mean,
                    stats.mapq_mean,
                    stats.edit_dist_mean,
                    stats.ani_mean,
                    stats.ani_std,
                    stats.ani_median,
                    <long long>stats.bases_covered,
                    <long long>stats.max_covered_bases,
                    stats.mean_covered_bases,
                    stats.mean_coverage,
                    stats.mean_coverage_trunc,
                    <long long>stats.mean_coverage_trunc_len,
                    stats.mean_coverage_covered,
                    <long long>stats.ref_length,
                    <long long>stats.bam_ref_length,
                    stats.breadth,
                    stats.exp_breadth,
                    stats.breadth_exp_ratio,
                    <long long>stats.n_bins,
                    stats.site_density,
                    stats.entropy,
                    stats.norm_entropy,
                    stats.gini,
                    stats.norm_gini,
                    stats.c_v,
                    stats.d_i,
                    stats.cov_evenness,
                    <long long>stats.tax_abund_read,
                    <long long>stats.tax_abund_aln,
                    <long long>stats.tax_abund_tad,
                    <long long>stats.n_reads_tad
                )
                refs_written += 1
                filtered_passed += 1

    if apply_filters:
        bf_nogil_logf_verbose(2, NULL, "[FILTER DEBUG] Total references with n_alns > 0: %d\n", filtered_total)
        bf_nogil_logf_verbose(2, NULL, "[FILTER DEBUG] References passing filters: %d\n", filtered_passed)
        bf_nogil_logf_verbose(2, NULL, "[FILTER DEBUG] Active filters:\n")

        if filters.enable_min_avg_read_ani:
            bf_nogil_logf_verbose(2, NULL, "  min_avg_read_ani: %.4f\n", filters.min_avg_read_ani)
        else:
            bf_nogil_logf_verbose(2, NULL, "  min_avg_read_ani: DISABLED\n")

        if filters.enable_min_expected_breadth_ratio:
            bf_nogil_logf_verbose(2, NULL, "  min_expected_breadth_ratio: %.4f\n", filters.min_expected_breadth_ratio)
        else:
            bf_nogil_logf_verbose(2, NULL, "  min_expected_breadth_ratio: DISABLED\n")

        if filters.enable_min_breadth:
            bf_nogil_logf_verbose(2, NULL, "  min_breadth: %.4f\n", filters.min_breadth)
        else:
            bf_nogil_logf_verbose(2, NULL, "  min_breadth: DISABLED\n")

        if filters.enable_min_coverage_evenness:
            bf_nogil_logf_verbose(2, NULL, "  min_coverage_evenness: %.4f\n", filters.min_coverage_evenness)
        else:
            bf_nogil_logf_verbose(2, NULL, "  min_coverage_evenness: DISABLED\n")

        if filters.enable_max_coeff_var:
            bf_nogil_logf_notime(NULL, "  max_coeff_var: %.4f\n", filters.max_coeff_var)
        else:
            bf_nogil_logf_notime(NULL, "  max_coeff_var: DISABLED\n")

        if filters.enable_min_coverage_mean:
            bf_nogil_logf_notime(NULL, "  min_coverage_mean: %.4f\n", filters.min_coverage_mean)
        else:
            bf_nogil_logf_notime(NULL, "  min_coverage_mean: DISABLED\n")

        if filters.enable_min_norm_entropy:
            bf_nogil_logf_notime(NULL, "  min_norm_entropy: %.4f\n", filters.min_norm_entropy)
        else:
            bf_nogil_logf_notime(NULL, "  min_norm_entropy: DISABLED\n")

        if filters.enable_max_norm_gini:
            bf_nogil_logf_notime(NULL, "  max_norm_gini: %.4f\n", filters.max_norm_gini)
        else:
            bf_nogil_logf_notime(NULL, "  max_norm_gini: DISABLED\n")

    return refs_written


cdef int write_stats_to_gzfile(
    gzFile gzfp,
    RefStats* global_ref_stats,
    sam_hdr_t* header,
    int n_refs,
    FilterConditions* filters,
    bint apply_filters
) except -1 nogil:
    cdef int i
    cdef RefStats* stats
    cdef int refs_written = 0
    cdef const char* ref_name

    cdef const char* OUTPUT_HEADER = (
        "reference\t"
        "n_reads\t"
        "n_alns\t"
        "read_length_mean\t"
        "read_length_std\t"
        "read_length_min\t"
        "read_length_max\t"
        "read_length_median\t"
        "read_length_mode\t"
        "gc_content_mean\t"
        "gc_content_std\t"
        "gc_content_total\t"
        "dust_mean\t"
        "dust_std\t"
        "read_aligned_length\t"
        "read_aln_score\t"
        "mapping_quality\t"
        "edit_distances\t"
        "read_ani_mean\t"
        "read_ani_std\t"
        "read_ani_median\t"
        "bases_covered\t"
        "max_covered_bases\t"
        "mean_covered_bases\t"
        "coverage_mean\t"
        "coverage_mean_trunc\t"
        "coverage_mean_trunc_len\t"
        "coverage_covered_mean\t"
        "reference_length\t"
        "bam_reference_length\t"
        "breadth\t"
        "exp_breadth\t"
        "breadth_exp_ratio\t"
        "n_bins\t"
        "site_density\t"
        "entropy\t"
        "norm_entropy\t"
        "gini\t"
        "norm_gini\t"
        "c_v\t"
        "d_i\t"
        "cov_evenness\t"
        "tax_abund_read\t"
        "tax_abund_aln\t"
        "tax_abund_tad\t"
        "n_reads_tad\n"
    )

    cdef const char* OUTPUT_FORMAT = (
        "%s\t"
        "%lld\t"
        "%lld\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%d\t"
        "%d\t"
        "%d\t"
        "%d\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%lld\t"
        "%lld\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%lld\t"
        "%.4f\t"
        "%lld\t"
        "%lld\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%lld\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%.4f\t"
        "%lld\t"
        "%lld\t"
        "%lld\t"
        "%lld\n"
    )

    gzprintf(gzfp, OUTPUT_HEADER)

    cdef int filtered_total = 0
    cdef int filtered_passed = 0

    for i in range(n_refs):
        stats = &global_ref_stats[i]
        if stats.n_alns > 0:
            filtered_total += 1
            ref_name = get_ref_name(header, i)
            if not apply_filters:
                gzprintf(
                    gzfp,
                    OUTPUT_FORMAT,
                    ref_name,
                    <long long>stats.n_reads,
                    <long long>stats.n_alns,
                    stats.read_length_mean,
                    stats.read_length_std,
                    stats.min_read_length,
                    stats.max_read_length,
                    stats.read_length_median,
                    stats.read_length_mode,
                    stats.gc_content_mean,
                    stats.gc_content_std,
                    stats.gc_content_total,
                    stats.dust_mean,
                    stats.dust_std,
                    stats.aligned_length_mean,
                    stats.aln_score_mean,
                    stats.mapq_mean,
                    stats.edit_dist_mean,
                    stats.ani_mean,
                    stats.ani_std,
                    stats.ani_median,
                    <long long>stats.bases_covered,
                    <long long>stats.max_covered_bases,
                    stats.mean_covered_bases,
                    stats.mean_coverage,
                    stats.mean_coverage_trunc,
                    <long long>stats.mean_coverage_trunc_len,
                    stats.mean_coverage_covered,
                    <long long>stats.ref_length,
                    <long long>stats.bam_ref_length,
                    stats.breadth,
                    stats.exp_breadth,
                    stats.breadth_exp_ratio,
                    <long long>stats.n_bins,
                    stats.site_density,
                    stats.entropy,
                    stats.norm_entropy,
                    stats.gini,
                    stats.norm_gini,
                    stats.c_v,
                    stats.d_i,
                    stats.cov_evenness,
                    <long long>stats.tax_abund_read,
                    <long long>stats.tax_abund_aln,
                    <long long>stats.tax_abund_tad,
                    <long long>stats.n_reads_tad
                )
                refs_written += 1
                filtered_passed += 1
            elif passes_filters(stats, filters):
                gzprintf(
                    gzfp,
                    OUTPUT_FORMAT,
                    ref_name,
                    <long long>stats.n_reads,
                    <long long>stats.n_alns,
                    stats.read_length_mean,
                    stats.read_length_std,
                    stats.min_read_length,
                    stats.max_read_length,
                    stats.read_length_median,
                    stats.read_length_mode,
                    stats.gc_content_mean,
                    stats.gc_content_std,
                    stats.gc_content_total,
                    stats.dust_mean,
                    stats.dust_std,
                    stats.aligned_length_mean,
                    stats.aln_score_mean,
                    stats.mapq_mean,
                    stats.edit_dist_mean,
                    stats.ani_mean,
                    stats.ani_std,
                    stats.ani_median,
                    <long long>stats.bases_covered,
                    <long long>stats.max_covered_bases,
                    stats.mean_covered_bases,
                    stats.mean_coverage,
                    stats.mean_coverage_trunc,
                    <long long>stats.mean_coverage_trunc_len,
                    stats.mean_coverage_covered,
                    <long long>stats.ref_length,
                    <long long>stats.bam_ref_length,
                    stats.breadth,
                    stats.exp_breadth,
                    stats.breadth_exp_ratio,
                    <long long>stats.n_bins,
                    stats.site_density,
                    stats.entropy,
                    stats.norm_entropy,
                    stats.gini,
                    stats.norm_gini,
                    stats.c_v,
                    stats.d_i,
                    stats.cov_evenness,
                    <long long>stats.tax_abund_read,
                    <long long>stats.tax_abund_aln,
                    <long long>stats.tax_abund_tad,
                    <long long>stats.n_reads_tad
                )
                refs_written += 1
                filtered_passed += 1

    return refs_written


cdef int write_output_files_complete(
    const char* output_c,
    const char* filtered_output_c,
    RefStats* global_ref_stats,
    sam_hdr_t* header,
    int n_refs,
    FilterConditions* filters
) except -1 nogil:
    cdef FILE* outfp = NULL
    cdef FILE* filtered_outfp = NULL
    cdef gzFile gzout = NULL
    cdef gzFile gzfilt = NULL
    cdef int refs_written = 0
    cdef int filtered_refs_written = 0

    # Write main output file (no filtering) - choose gz when requested
    if output_c != NULL:
        if c_str_endswith_gz(output_c):
            gzout = gzopen(output_c, "wb")
            if gzout != NULL:
                refs_written = write_stats_to_gzfile(gzout, global_ref_stats, header, n_refs, filters, False)
                gzclose(gzout)
        else:
            outfp = fopen(output_c, "w")
            if outfp != NULL:
                refs_written = write_stats_to_file(outfp, global_ref_stats, header, n_refs, filters, False)
                fclose(outfp)

    # Write filtered output file (with filtering)
    if filtered_output_c != NULL:
        if c_str_endswith_gz(filtered_output_c):
            gzfilt = gzopen(filtered_output_c, "wb")
            if gzfilt != NULL:
                filtered_refs_written = write_stats_to_gzfile(gzfilt, global_ref_stats, header, n_refs, filters, True)
                gzclose(gzfilt)
        else:
            filtered_outfp = fopen(filtered_output_c, "w")
            if filtered_outfp != NULL:
                filtered_refs_written = write_stats_to_file(filtered_outfp, global_ref_stats, header, n_refs, filters, True)
                fclose(filtered_outfp)

    return 0
