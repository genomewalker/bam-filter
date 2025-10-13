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

"""Processing statistics tracking and reporting.

Tracks comprehensive statistics through all processing stages: initial BAM reading,
quality filtering, EM algorithm, probability filtering, graph analysis, and final output.
"""

from libc.stdint cimport int64_t, int32_t
from libc.string cimport memset
from .processor cimport ProcessingStats, MemoryPool

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil


cdef void init_processing_stats(ProcessingStats* stats) noexcept nogil:
    """Initialize all statistics to zero."""
    memset(stats, 0, sizeof(ProcessingStats))


cdef void update_initial_stats(ProcessingStats* stats, int64_t alignments, 
                               int64_t reads, int64_t references) noexcept nogil:
    """Update initial BAM reading statistics."""
    stats.initial_total_alignments = alignments
    stats.initial_total_reads = reads
    stats.initial_total_references = references


cdef void update_quality_filter_stats(ProcessingStats* stats, int64_t alignments,
                                      int64_t reads, int64_t references) noexcept nogil:
    """Update quality filtering (ANI/length) statistics."""
    stats.post_quality_alignments = alignments
    stats.post_quality_reads = reads
    stats.post_quality_references = references
    stats.filtered_quality_alignments = stats.initial_total_alignments - alignments


cdef void update_em_stats(ProcessingStats* stats, int32_t iterations,
                         bint converged, double likelihood) noexcept nogil:
    """Update EM algorithm statistics."""
    stats.em_iterations = iterations
    stats.em_converged = converged
    stats.em_final_likelihood = likelihood


cdef void update_probability_filter_stats(ProcessingStats* stats, int64_t alignments,
                                          int64_t reads, int64_t references) noexcept nogil:
    """Update probability filtering statistics."""
    stats.post_probability_alignments = alignments
    stats.post_probability_reads = reads
    stats.post_probability_references = references
    stats.filtered_probability_alignments = stats.post_quality_alignments - alignments


cdef void update_graph_stats(ProcessingStats* stats, int64_t references,
                             int64_t patterns) noexcept nogil:
    """Update graph analysis statistics."""
    stats.graph_analysis_references = references
    stats.graph_patterns_computed = patterns


cdef void update_unified_filter_stats(ProcessingStats* stats,
                                      int64_t alignments, int64_t reads, int64_t references,
                                      int64_t coverage_only, int64_t info_only, int64_t both,
                                      int64_t align_cov, int64_t align_info, int64_t align_both) noexcept nogil:
    """Update unified filtering statistics."""
    stats.post_unified_alignments = alignments
    stats.post_unified_reads = reads
    stats.post_unified_references = references
    
    stats.filtered_coverage_only = coverage_only
    stats.filtered_information_only = info_only
    stats.filtered_both_criteria = both
    
    stats.alignments_removed_coverage = align_cov
    stats.alignments_removed_information = align_info
    stats.alignments_removed_both = align_both


cdef void update_final_output_stats(ProcessingStats* stats, int64_t alignments,
                                    int64_t reads, int64_t references) noexcept nogil:
    """Update final output statistics."""
    stats.final_alignments_written = alignments
    stats.final_reads_written = reads
    stats.final_references_written = references


cdef void calculate_summary_metrics(ProcessingStats* stats) noexcept nogil:
    """Calculate summary retention percentages."""
    if stats.initial_total_alignments > 0:
        stats.overall_alignment_retention = (100.0 * <double>stats.final_alignments_written / 
                                            <double>stats.initial_total_alignments)
    else:
        stats.overall_alignment_retention = 0.0
    
    if stats.initial_total_reads > 0:
        stats.overall_read_retention = (100.0 * <double>stats.final_reads_written / 
                                       <double>stats.initial_total_reads)
    else:
        stats.overall_read_retention = 0.0
    
    if stats.initial_total_references > 0:
        stats.overall_reference_retention = (100.0 * <double>stats.final_references_written / 
                                            <double>stats.initial_total_references)
    else:
        stats.overall_reference_retention = 0.0


cdef void print_processing_stats(ProcessingStats* stats) noexcept nogil:
    """Print comprehensive formatted statistics report."""
    cdef const char* tag = b"STATS"

    bf_nogil_logf_notime(tag, "")
    bf_nogil_logf_notime(tag, "================================================================================\n")
    bf_nogil_logf_notime(tag, "                    COMPREHENSIVE PROCESSING STATISTICS\n")
    bf_nogil_logf_notime(tag, "================================================================================\n")

    bf_nogil_logf_notime(tag, "[STAGE 1] Initial BAM Reading\n")
    bf_nogil_logf_notime(tag, "  Total alignments:     %12lld\n", <long long>stats.initial_total_alignments)
    bf_nogil_logf_notime(tag, "  Unique reads:         %12lld\n", <long long>stats.initial_total_reads)
    bf_nogil_logf_notime(tag, "  References with data: %12lld\n\n", <long long>stats.initial_total_references)

    if stats.post_quality_alignments > 0:
        bf_nogil_logf_notime(tag, "[STAGE 2] Quality Filtering (ANI >= %.1f%%, Length filters)\n", 90.0)
        bf_nogil_logf_notime(
            tag,
            "  Alignments kept:      %12lld (%.1f%%)\n",
            <long long>stats.post_quality_alignments,
            100.0 * <double>stats.post_quality_alignments / <double>stats.initial_total_alignments,
        )
        bf_nogil_logf_notime(tag, "  Alignments filtered:  %12lld\n", <long long>stats.filtered_quality_alignments)
        bf_nogil_logf_notime(tag, "  Unique reads kept:    %12lld\n", <long long>stats.post_quality_reads)
        bf_nogil_logf_notime(tag, "  References kept:      %12lld\n\n", <long long>stats.post_quality_references)

    bf_nogil_logf_notime(tag, "[STAGE 3] EM Algorithm (Read Reassignment)\n")
    bf_nogil_logf_notime(tag, "  Iterations:           %12d\n", stats.em_iterations)
    bf_nogil_logf_notime(tag, "  Converged:            %12s\n", b"YES" if stats.em_converged else b"NO")
    bf_nogil_logf_notime(tag, "  Final log-likelihood: %12.6f\n\n", stats.em_final_likelihood)

    if stats.post_probability_alignments > 0:
        bf_nogil_logf_notime(tag, "[STAGE 4] Probability Filtering\n")
        if stats.post_quality_alignments > 0:
            bf_nogil_logf_notime(
                tag,
                "  Alignments kept:      %12lld (%.1f%%)\n",
                <long long>stats.post_probability_alignments,
                100.0 * <double>stats.post_probability_alignments / <double>stats.post_quality_alignments,
            )
        else:
            bf_nogil_logf_notime(
                tag,
                "  Alignments kept:      %12lld\n",
                <long long>stats.post_probability_alignments,
            )
        bf_nogil_logf_notime(tag, "  Alignments filtered:  %12lld\n", <long long>stats.filtered_probability_alignments)
        bf_nogil_logf_notime(tag, "  Unique reads:         %12lld\n", <long long>stats.post_probability_reads)
        bf_nogil_logf_notime(tag, "  References:           %12lld\n\n", <long long>stats.post_probability_references)

    if stats.graph_analysis_references >= 0:
        bf_nogil_logf_notime(tag, "[STAGE 5] Graph Analysis (Network Patterns)\n")
        bf_nogil_logf_notime(tag, "  References analyzed:  %12lld\n", <long long>stats.graph_analysis_references)
        if stats.graph_analysis_references > 0:
            bf_nogil_logf_notime(
                tag,
                "  Patterns computed:    %12lld (%.1f%%)\n\n",
                <long long>stats.graph_patterns_computed,
                100.0 * <double>stats.graph_patterns_computed / <double>stats.graph_analysis_references,
            )
        else:
            bf_nogil_logf_notime(
                tag,
                "  Patterns computed:    %12lld\n\n",
                <long long>stats.graph_patterns_computed,
            )
    else:
        bf_nogil_logf_notime(tag, "[STAGE 5] Graph Analysis (Network Patterns)\n")
        bf_nogil_logf_notime(tag, "  Status:              skipped (graph export disabled)\n\n")

    if stats.filtered_coverage_only >= 0:
        bf_nogil_logf_notime(tag, "[STAGE 6] Unified Filtering (Coverage + Information)\n")
        if stats.post_probability_references > 0:
            bf_nogil_logf_notime(
                tag,
                "  References kept:      %12lld (%.1f%%)\n",
                <long long>stats.post_unified_references,
                100.0 * <double>stats.post_unified_references / <double>stats.post_probability_references,
            )
        else:
            bf_nogil_logf_notime(tag, "  References kept:      %12lld\n", <long long>stats.post_unified_references)
        bf_nogil_logf_notime(tag, "  References filtered:\n")
        bf_nogil_logf_notime(tag, "    Coverage only:      %12lld\n", <long long>stats.filtered_coverage_only)
        bf_nogil_logf_notime(tag, "    Information only:   %12lld\n", <long long>stats.filtered_information_only)
        bf_nogil_logf_notime(tag, "    Both criteria:      %12lld\n\n", <long long>stats.filtered_both_criteria)

        if stats.post_probability_alignments > 0:
            bf_nogil_logf_notime(
                tag,
                "  Alignments kept:      %12lld (%.1f%%)\n",
                <long long>stats.post_unified_alignments,
                100.0 * <double>stats.post_unified_alignments / <double>stats.post_probability_alignments,
            )
        else:
            bf_nogil_logf_notime(tag, "  Alignments kept:      %12lld\n", <long long>stats.post_unified_alignments)
        bf_nogil_logf_notime(tag, "  Alignments removed:\n")
        bf_nogil_logf_notime(tag, "    From coverage:      %12lld\n", <long long>stats.alignments_removed_coverage)
        bf_nogil_logf_notime(tag, "    From information:   %12lld\n", <long long>stats.alignments_removed_information)
        bf_nogil_logf_notime(tag, "    From both:          %12lld\n\n", <long long>stats.alignments_removed_both)

        if stats.post_probability_reads > 0:
            bf_nogil_logf_notime(
                tag,
                "  Unique reads kept:    %12lld (%.1f%%)\n\n",
                <long long>stats.post_unified_reads,
                100.0 * <double>stats.post_unified_reads / <double>stats.post_probability_reads,
            )
        else:
            bf_nogil_logf_notime(tag, "  Unique reads kept:    %12lld\n\n", <long long>stats.post_unified_reads)
    else:
        bf_nogil_logf_notime(tag, "[STAGE 6] Unified Filtering (Coverage + Information)\n")
        bf_nogil_logf_notime(tag, "  Status:              skipped (cluster-aware filtering disabled)\n")
        bf_nogil_logf_notime(
            tag,
            "  References forwarded: %12lld\n",
            <long long>stats.post_unified_references,
        )
        bf_nogil_logf_notime(
            tag,
            "  Alignments forwarded: %12lld\n",
            <long long>stats.post_unified_alignments,
        )
        bf_nogil_logf_notime(
            tag,
            "  Unique reads forwarded: %12lld\n\n",
                <long long>stats.post_unified_reads,
        )

    bf_nogil_logf_notime(tag, "[STAGE 7] Final Output BAM\n")
    bf_nogil_logf_notime(tag, "  Alignments written:   %12lld\n", <long long>stats.final_alignments_written)
    bf_nogil_logf_notime(tag, "  Unique reads written: %12lld\n", <long long>stats.final_reads_written)
    bf_nogil_logf_notime(tag, "  References written:   %12lld\n\n", <long long>stats.final_references_written)

    bf_nogil_logf_notime(tag, "================================================================================\n")
    bf_nogil_logf_notime(tag, "                           OVERALL RETENTION\n")
    bf_nogil_logf_notime(tag, "================================================================================\n")
    bf_nogil_logf_notime(
        tag,
        "  Alignments:  %12lld / %12lld (%.2f%%)\n",
        <long long>stats.final_alignments_written,
        <long long>stats.initial_total_alignments,
        stats.overall_alignment_retention,
    )
    bf_nogil_logf_notime(
        tag,
        "  Reads:       %12lld / %12lld (%.2f%%)\n",
        <long long>stats.final_reads_written,
        <long long>stats.initial_total_reads,
        stats.overall_read_retention,
    )
    bf_nogil_logf_notime(
        tag,
        "  References:  %12lld / %12lld (%.2f%%)\n",
        <long long>stats.final_references_written,
        <long long>stats.initial_total_references,
        stats.overall_reference_retention,
    )
    bf_nogil_logf_notime(tag, "================================================================================\n\n")
