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

"""Generic column-based filtering system for BAM statistics.

This module provides a flexible filtering system that can filter on any
of the 45 exported TSV columns.
"""

from libc.stdlib cimport malloc, free, realloc
from libc.math cimport INFINITY, isnan
from bam_filter.stats cimport RefStats

# Maximum number of filters we expect (can be increased)
cdef int MAX_FILTERS = 50

cdef struct ColumnFilter:
    int column_index
    double min_value
    double max_value
    bint is_active

cdef struct GenericFilters:
    ColumnFilter* filters
    int n_filters
    int capacity


cdef GenericFilters* create_generic_filters(int initial_capacity) noexcept nogil:
    """Create a new GenericFilters structure.

    Args:
        initial_capacity: Initial number of filter slots to allocate

    Returns:
        Pointer to allocated GenericFilters, or NULL on failure
    """
    cdef GenericFilters* gf = <GenericFilters*>malloc(sizeof(GenericFilters))
    if gf == NULL:
        return NULL

    gf.capacity = initial_capacity if initial_capacity > 0 else 10
    gf.n_filters = 0
    gf.filters = <ColumnFilter*>malloc(gf.capacity * sizeof(ColumnFilter))

    if gf.filters == NULL:
        free(gf)
        return NULL

    return gf


cdef void destroy_generic_filters(GenericFilters* gf) noexcept nogil:
    """Free memory used by GenericFilters."""
    if gf != NULL:
        if gf.filters != NULL:
            free(gf.filters)
        free(gf)


cdef int add_filter(GenericFilters* gf, int column_index, double min_val, double max_val) noexcept nogil:
    """Add a filter rule to the GenericFilters.

    Args:
        gf: GenericFilters structure
        column_index: Column index (0-44)
        min_val: Minimum value (-INFINITY for no lower bound)
        max_val: Maximum value (+INFINITY for no upper bound)

    Returns:
        0 on success, -1 on failure
    """
    cdef int new_capacity
    cdef ColumnFilter* new_filters

    if gf == NULL:
        return -1

    # Resize if needed
    if gf.n_filters >= gf.capacity:
        new_capacity = gf.capacity * 2
        new_filters = <ColumnFilter*>realloc(
            gf.filters, new_capacity * sizeof(ColumnFilter)
        )
        if new_filters == NULL:
            return -1
        gf.filters = new_filters
        gf.capacity = new_capacity

    # Add the filter
    gf.filters[gf.n_filters].column_index = column_index
    gf.filters[gf.n_filters].min_value = min_val
    gf.filters[gf.n_filters].max_value = max_val
    gf.filters[gf.n_filters].is_active = True
    gf.n_filters += 1

    return 0


cdef double get_column_value(RefStats* stats, int column_index) noexcept nogil:
    """Extract the value for a specific column from RefStats.

    Column indices correspond to TSV output order:
    0: reference_name (not numeric, returns 0)
    1: read_count (n_reads)
    2: alignment_count (n_alns)
    ... and so on

    Args:
        stats: RefStats structure
        column_index: Column index (0-44)

    Returns:
        Column value as double, or 0 for non-numeric columns
    """
    # Column mapping based on TSV output order
    if column_index == 0:  # reference_name - not numeric
        return 0.0
    elif column_index == 1:  # read_count
        return <double>stats.n_reads
    elif column_index == 2:  # alignment_count
        return <double>stats.n_alns
    elif column_index == 3:  # read_length_mean
        return stats.read_length_mean
    elif column_index == 4:  # read_length_std
        return stats.read_length_std
    elif column_index == 5:  # read_length_min
        return <double>stats.min_read_length
    elif column_index == 6:  # read_length_max
        return <double>stats.max_read_length
    elif column_index == 7:  # read_length_median
        return <double>stats.read_length_median
    elif column_index == 8:  # read_length_mode
        return <double>stats.read_length_mode
    elif column_index == 9:  # gc_content_mean (backward compat: maps to read_gc)
        return stats.read_gc_content_mean
    elif column_index == 10:  # gc_content_std (backward compat: maps to read_gc)
        return stats.read_gc_content_std
    elif column_index == 11:  # gc_content_total (backward compat: maps to read_gc)
        return stats.read_gc_content_total
    elif column_index == 12:  # dust_mean
        return stats.dust_mean
    elif column_index == 13:  # dust_std
        return stats.dust_std
    elif column_index == 14:  # read_aligned_length_mean
        return stats.aligned_length_mean
    elif column_index == 15:  # read_alignment_score_mean
        return stats.aln_score_mean
    elif column_index == 16:  # mapping_quality_mean
        return stats.mapq_mean
    elif column_index == 17:  # edit_distance_mean
        return stats.edit_dist_mean
    elif column_index == 18:  # read_ani_mean
        return stats.ani_mean
    elif column_index == 19:  # read_ani_std
        return stats.ani_std
    elif column_index == 20:  # read_ani_median
        return stats.ani_median
    elif column_index == 21:  # bases_covered
        return <double>stats.bases_covered
    elif column_index == 22:  # bases_covered_max
        return <double>stats.max_covered_bases
    elif column_index == 23:  # bases_covered_mean
        return stats.mean_covered_bases
    elif column_index == 24:  # coverage_mean
        return stats.mean_coverage
    elif column_index == 25:  # coverage_mean_trimmed
        return stats.mean_coverage_trunc
    elif column_index == 26:  # coverage_mean_trimmed_length
        return <double>stats.mean_coverage_trunc_len
    elif column_index == 27:  # coverage_mean_covered_only
        return stats.mean_coverage_covered
    elif column_index == 28:  # reference_length
        return <double>stats.ref_length
    elif column_index == 29:  # reference_length_bam
        return <double>stats.bam_ref_length
    elif column_index == 30:  # breadth
        return stats.breadth
    elif column_index == 31:  # breadth_expected
        return stats.exp_breadth
    elif column_index == 32:  # breadth_expected_ratio
        return stats.breadth_exp_ratio
    elif column_index == 33:  # bin_count
        return <double>stats.n_bins
    elif column_index == 34:  # site_density
        return stats.site_density
    elif column_index == 35:  # spatial_entropy
        return stats.spatial_entropy
    elif column_index == 36:  # spatial_entropy_normalized
        return stats.norm_spatial_entropy
    elif column_index == 37:  # gini_coefficient
        return stats.gini
    elif column_index == 38:  # gini_coefficient_normalized
        return stats.norm_gini
    elif column_index == 39:  # coefficient_of_variation
        return stats.c_v
    elif column_index == 40:  # diversity_index
        return stats.d_i
    elif column_index == 41:  # coverage_evenness
        return stats.cov_evenness
    elif column_index == 42:  # abundance_read_based
        return stats.tax_abund_read
    elif column_index == 43:  # abundance_alignment_based
        return stats.tax_abund_aln
    elif column_index == 44:  # abundance_tad
        return stats.tax_abund_tad
    elif column_index == 45:  # read_count_tad
        return <double>stats.n_reads_tad
    else:
        return 0.0


cdef bint passes_generic_filters(RefStats* stats, GenericFilters* gf) noexcept nogil:
    """Check if a reference passes all filter conditions.

    Args:
        stats: RefStats to check
        gf: GenericFilters with filter rules

    Returns:
        True if passes all filters, False otherwise
    """
    if gf == NULL or gf.n_filters == 0:
        return True  # No filters means pass

    cdef int i
    cdef double value
    cdef ColumnFilter* f

    for i in range(gf.n_filters):
        f = &gf.filters[i]
        if not f.is_active:
            continue

        value = get_column_value(stats, f.column_index)

        # Skip NaN values (treat as passing)
        if isnan(value):
            continue

        # Check bounds
        if value < f.min_value or value > f.max_value:
            return False

    return True
