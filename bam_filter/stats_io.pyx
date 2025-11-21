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
Statistics output writer using unified TSV writer.

This module provides a clean interface for writing BAM filter statistics
to TSV files with optional compression.
"""

from libc.string cimport strlen
from bam_filter.stats cimport RefStats, FilterConditions
from bam_filter.generic_filters cimport GenericFilters, passes_generic_filters
from bam_filter.processor cimport sam_hdr_t, sam_hdr_tid2name
from bam_filter.tsv_writer cimport (
    TSVWriter, tsv_writer_open, tsv_writer_close, tsv_writer_flush,
    tsv_row_start, tsv_append_string, tsv_append_int, tsv_append_float, tsv_row_end
)

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_verbose(int level, const char* tag, const char* fmt, ...) nogil

cdef const char* UNKNOWN_REF_NAME = "<unknown>"

cdef inline const char* get_ref_name(sam_hdr_t* header, int tid) noexcept nogil:
    """Get reference name from BAM header, or return <unknown> if not found."""
    cdef const char* name = sam_hdr_tid2name(header, tid)
    if name == NULL:
        return UNKNOWN_REF_NAME
    return name


cdef int write_stats_header(TSVWriter* writer) except -1 nogil:
    """Write the TSV header row with all column names."""
    tsv_row_start(writer)
    tsv_append_string(writer, "reference")
    tsv_append_string(writer, "n_reads")
    tsv_append_string(writer, "n_alns")
    tsv_append_string(writer, "read_length_mean")
    tsv_append_string(writer, "read_length_std")
    tsv_append_string(writer, "read_length_min")
    tsv_append_string(writer, "read_length_max")
    tsv_append_string(writer, "read_length_median")
    tsv_append_string(writer, "read_length_mode")
    tsv_append_string(writer, "gc_content_mean")
    tsv_append_string(writer, "gc_content_std")
    tsv_append_string(writer, "gc_content_total")
    tsv_append_string(writer, "dust_mean")
    tsv_append_string(writer, "dust_std")
    tsv_append_string(writer, "read_aligned_length")
    tsv_append_string(writer, "read_aln_score")
    tsv_append_string(writer, "mapping_quality")
    tsv_append_string(writer, "edit_distances")
    tsv_append_string(writer, "read_ani_mean")
    tsv_append_string(writer, "read_ani_std")
    tsv_append_string(writer, "read_ani_median")
    tsv_append_string(writer, "bases_covered")
    tsv_append_string(writer, "max_covered_bases")
    tsv_append_string(writer, "mean_covered_bases")
    tsv_append_string(writer, "coverage_mean")
    tsv_append_string(writer, "coverage_mean_trunc")
    tsv_append_string(writer, "coverage_mean_trunc_len")
    tsv_append_string(writer, "coverage_covered_mean")
    tsv_append_string(writer, "reference_length")
    tsv_append_string(writer, "bam_reference_length")
    tsv_append_string(writer, "breadth")
    tsv_append_string(writer, "exp_breadth")
    tsv_append_string(writer, "breadth_exp_ratio")
    tsv_append_string(writer, "n_bins")
    tsv_append_string(writer, "site_density")
    tsv_append_string(writer, "spatial_entropy")
    tsv_append_string(writer, "norm_spatial_entropy")
    tsv_append_string(writer, "gini")
    tsv_append_string(writer, "norm_gini")
    tsv_append_string(writer, "c_v")
    tsv_append_string(writer, "d_i")
    tsv_append_string(writer, "cov_evenness")
    tsv_append_string(writer, "tax_abund_read")
    tsv_append_string(writer, "tax_abund_aln")
    tsv_append_string(writer, "tax_abund_tad")
    tsv_append_string(writer, "n_reads_tad")
    tsv_row_end(writer)
    return 0


cdef int write_stats_row(TSVWriter* writer, RefStats* stats, const char* ref_name) except -1 nogil:
    """Write a single statistics row to the TSV file."""
    tsv_row_start(writer)
    tsv_append_string(writer, ref_name)
    tsv_append_int(writer, stats.n_reads)
    tsv_append_int(writer, stats.n_alns)
    tsv_append_float(writer, stats.read_length_mean, 4)
    tsv_append_float(writer, stats.read_length_std, 4)
    tsv_append_int(writer, stats.min_read_length)
    tsv_append_int(writer, stats.max_read_length)
    tsv_append_int(writer, stats.read_length_median)
    tsv_append_int(writer, stats.read_length_mode)
    tsv_append_float(writer, stats.gc_content_mean, 4)
    tsv_append_float(writer, stats.gc_content_std, 4)
    tsv_append_float(writer, stats.gc_content_total, 4)
    tsv_append_float(writer, stats.dust_mean, 4)
    tsv_append_float(writer, stats.dust_std, 4)
    tsv_append_float(writer, stats.aligned_length_mean, 4)
    tsv_append_float(writer, stats.aln_score_mean, 4)
    tsv_append_float(writer, stats.mapq_mean, 4)
    tsv_append_float(writer, stats.edit_dist_mean, 4)
    tsv_append_float(writer, stats.ani_mean, 4)
    tsv_append_float(writer, stats.ani_std, 4)
    tsv_append_float(writer, stats.ani_median, 4)
    tsv_append_int(writer, stats.bases_covered)
    tsv_append_int(writer, stats.max_covered_bases)
    tsv_append_float(writer, stats.mean_covered_bases, 4)
    tsv_append_float(writer, stats.mean_coverage, 4)
    tsv_append_float(writer, stats.mean_coverage_trunc, 4)
    tsv_append_int(writer, stats.mean_coverage_trunc_len)
    tsv_append_float(writer, stats.mean_coverage_covered, 4)
    tsv_append_int(writer, stats.ref_length)
    tsv_append_int(writer, stats.bam_ref_length)
    tsv_append_float(writer, stats.breadth, 4)
    tsv_append_float(writer, stats.exp_breadth, 4)
    tsv_append_float(writer, stats.breadth_exp_ratio, 4)
    tsv_append_int(writer, stats.n_bins)
    tsv_append_float(writer, stats.site_density, 4)
    tsv_append_float(writer, stats.spatial_entropy, 4)
    tsv_append_float(writer, stats.norm_spatial_entropy, 4)
    tsv_append_float(writer, stats.gini, 4)
    tsv_append_float(writer, stats.norm_gini, 4)
    tsv_append_float(writer, stats.c_v, 4)
    tsv_append_float(writer, stats.d_i, 4)
    tsv_append_float(writer, stats.cov_evenness, 4)
    tsv_append_int(writer, stats.tax_abund_read)
    tsv_append_int(writer, stats.tax_abund_aln)
    tsv_append_int(writer, stats.tax_abund_tad)
    tsv_append_int(writer, stats.n_reads_tad)
    tsv_row_end(writer)
    return 0


cdef int write_stats_to_writer(
    TSVWriter* writer,
    RefStats* global_ref_stats,
    sam_hdr_t* header,
    int n_refs,
    GenericFilters* gfilters,
    bint apply_filters
) except -1 nogil:
    """
    Write statistics to a TSV writer with optional filtering.

    Parameters
    ----------
    writer : TSVWriter*
        Opened TSV writer
    global_ref_stats : RefStats*
        Array of reference statistics
    header : sam_hdr_t*
        BAM header for reference names
    n_refs : int
        Number of references
    gfilters : GenericFilters*
        Generic filters to apply (can be NULL)
    apply_filters : bint
        Whether to apply filters

    Returns
    -------
    int
        Number of references written, or -1 on error
    """
    cdef int i
    cdef RefStats* stats
    cdef int refs_written = 0
    cdef const char* ref_name
    cdef int filtered_total = 0
    cdef int filtered_passed = 0

    # Write header
    write_stats_header(writer)

    # Write data rows
    for i in range(n_refs):
        stats = &global_ref_stats[i]
        if stats.n_alns > 0:
            filtered_total += 1
            ref_name = get_ref_name(header, i)

            if not apply_filters:
                write_stats_row(writer, stats, ref_name)
                refs_written += 1
                filtered_passed += 1
            elif passes_generic_filters(stats, gfilters):
                write_stats_row(writer, stats, ref_name)
                refs_written += 1
                filtered_passed += 1

    # Log filtering statistics
    if apply_filters:
        bf_nogil_logf_verbose(2, NULL, "[FILTER] Total references with n_alns > 0: %d\n", filtered_total)
        bf_nogil_logf_verbose(2, NULL, "[FILTER] References passing filters: %d\n", filtered_passed)
        if gfilters != NULL and gfilters.n_filters > 0:
            bf_nogil_logf_verbose(2, NULL, "[FILTER] Active filters: %d\n", gfilters.n_filters)

    return refs_written


cdef int write_output_files_complete(
    const char* output_c,
    const char* filtered_output_c,
    RefStats* global_ref_stats,
    sam_hdr_t* header,
    int n_refs,
    FilterConditions* filters,
    GenericFilters* gfilters
) except -1 nogil:
    """
    Write output files with default compression settings.
    This is a convenience wrapper that uses fast compression (level 1, 8 threads).
    """
    return write_output_files_complete_fast(
        output_c, filtered_output_c, global_ref_stats,
        header, n_refs, filters, gfilters, 1, 8
    )


cdef int write_output_files_complete_fast(
    const char* output_c,
    const char* filtered_output_c,
    RefStats* global_ref_stats,
    sam_hdr_t* header,
    int n_refs,
    FilterConditions* filters,
    GenericFilters* gfilters,
    int compression_level,
    int compression_threads
) except -1 nogil:
    """
    Write statistics output files with configurable compression.

    This is the main entry point for writing statistics. It handles:
    - Automatic compression detection (based on .gz extension)
    - Pigz multi-threaded compression (falls back to zlib if unavailable)
    - Filtered and unfiltered output
    - Ultra-fast TSV writing (no printf overhead)

    Parameters
    ----------
    output_c : const char*
        Main output file path (can be NULL to skip)
    filtered_output_c : const char*
        Filtered output file path (can be NULL to skip)
    global_ref_stats : RefStats*
        Array of reference statistics
    header : sam_hdr_t*
        BAM header for reference names
    n_refs : int
        Number of references
    filters : FilterConditions*
        Legacy filter conditions (deprecated, use gfilters)
    gfilters : GenericFilters*
        Generic filters to apply
    compression_level : int
        Compression level (1=fastest, 9=best). Recommend 1 for speed, 6 for size.
    compression_threads : int
        Number of threads for pigz compression

    Returns
    -------
    int
        0 on success, -1 on error
    """
    cdef TSVWriter* writer = NULL
    cdef int refs_written = 0

    # Write main output file (no filtering)
    if output_c != NULL:
        writer = tsv_writer_open(output_c, compression_level, compression_threads)
        if writer == NULL:
            return -1

        refs_written = write_stats_to_writer(writer, global_ref_stats, header, n_refs, gfilters, False)
        tsv_writer_close(writer)

        if refs_written < 0:
            return -1

    # Write filtered output file (with filtering)
    if filtered_output_c != NULL:
        writer = tsv_writer_open(filtered_output_c, compression_level, compression_threads)
        if writer == NULL:
            return -1

        refs_written = write_stats_to_writer(writer, global_ref_stats, header, n_refs, gfilters, True)
        tsv_writer_close(writer)

        if refs_written < 0:
            return -1

    return 0


# Legacy compatibility functions (deprecated but kept for API compatibility)
# These are maintained for backwards compatibility but should not be used in new code

from libc.stdio cimport FILE, fprintf
from libc.stdlib cimport malloc, free
from bam_filter.tsv_writer cimport COMPRESSION_NONE, COMPRESSION_ZLIB

cdef extern from "zlib.h":
    ctypedef void* gzFile
    int gzprintf(gzFile file, const char* format, ...) nogil

cdef int write_stats_to_file(
    FILE* fp,
    RefStats* global_ref_stats,
    sam_hdr_t* header,
    int n_refs,
    FilterConditions* filters,
    GenericFilters* gfilters,
    bint apply_filters
) except -1 nogil:
    """
    DEPRECATED: Legacy function for writing to FILE*.
    Kept for backwards compatibility. New code should use write_stats_to_writer.
    """
    # Create a temporary TSV writer that wraps the FILE*
    # This is inefficient but maintains compatibility
    cdef TSVWriter* writer = <TSVWriter*>malloc(sizeof(TSVWriter))
    if writer == NULL:
        return -1

    writer.file_handle = fp
    writer.gz_handle = NULL
    writer.compression_type = COMPRESSION_NONE
    writer.line_capacity = 8192
    writer.output_capacity = 67108864
    writer.line_pos = 0
    writer.output_pos = 0

    writer.line_buffer = <char*>malloc(writer.line_capacity)
    writer.output_buffer = <char*>malloc(writer.output_capacity)

    if writer.line_buffer == NULL or writer.output_buffer == NULL:
        if writer.line_buffer != NULL:
            free(writer.line_buffer)
        if writer.output_buffer != NULL:
            free(writer.output_buffer)
        free(writer)
        return -1

    cdef int result = write_stats_to_writer(writer, global_ref_stats, header, n_refs, gfilters, apply_filters)

    # Flush but don't close the file (caller owns it)
    tsv_writer_flush(writer)
    free(writer.line_buffer)
    free(writer.output_buffer)
    free(writer)

    return result


cdef int write_stats_to_gzfile(
    gzFile gzfp,
    RefStats* global_ref_stats,
    sam_hdr_t* header,
    int n_refs,
    FilterConditions* filters,
    GenericFilters* gfilters,
    bint apply_filters
) except -1 nogil:
    """
    DEPRECATED: Legacy function for writing to gzFile.
    Kept for backwards compatibility. New code should use write_stats_to_writer.
    """
    # Create a temporary TSV writer that wraps the gzFile
    cdef TSVWriter* writer = <TSVWriter*>malloc(sizeof(TSVWriter))
    if writer == NULL:
        return -1

    writer.file_handle = NULL
    writer.gz_handle = gzfp
    writer.compression_type = COMPRESSION_ZLIB
    writer.line_capacity = 8192
    writer.output_capacity = 67108864
    writer.line_pos = 0
    writer.output_pos = 0

    writer.line_buffer = <char*>malloc(writer.line_capacity)
    writer.output_buffer = <char*>malloc(writer.output_capacity)

    if writer.line_buffer == NULL or writer.output_buffer == NULL:
        if writer.line_buffer != NULL:
            free(writer.line_buffer)
        if writer.output_buffer != NULL:
            free(writer.output_buffer)
        free(writer)
        return -1

    cdef int result = write_stats_to_writer(writer, global_ref_stats, header, n_refs, gfilters, apply_filters)

    # Flush but don't close the file (caller owns it)
    tsv_writer_flush(writer)
    free(writer.line_buffer)
    free(writer.output_buffer)
    free(writer)

    return result
