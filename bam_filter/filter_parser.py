"""Filter specification parser for generic column-based filtering.

This module provides functionality to parse and validate filter specifications
for BAM filtering operations. Filters use the format: column:min:max

Examples:
    read_ani_mean:90:       # >= 90
    breadth::0.95           # <= 0.95
    coverage_mean:5:1000    # between 5 and 1000
"""

from typing import Dict, Optional, Tuple, Any
from enum import Enum


class ColumnType(Enum):
    """Data types for filterable columns."""
    INT = "int"
    FLOAT = "float"
    STR = "str"


class ColumnSpec:
    """Specification for a filterable column including name, type, and valid range."""

    def __init__(
        self,
        old_name: str,
        new_name: str,
        dtype: ColumnType,
        min_value: Optional[float] = None,
        max_value: Optional[float] = None,
        filterable: bool = True
    ):
        self.old_name = old_name
        self.new_name = new_name
        self.dtype = dtype
        self.min_value = min_value
        self.max_value = max_value
        self.filterable = filterable

    def __repr__(self):
        bounds = []
        if self.min_value is not None:
            bounds.append(f"min={self.min_value}")
        if self.max_value is not None:
            bounds.append(f"max={self.max_value}")
        bounds_str = f" ({', '.join(bounds)})" if bounds else ""
        return f"{self.new_name} [{self.dtype.value}]{bounds_str}"


# Column specifications: old_name -> (new_name, dtype, min, max, filterable)
COLUMN_SPECS = {
    # Reference & Basic Counts
    'reference': ColumnSpec('reference', 'reference_name', ColumnType.STR, filterable=False),
    'n_reads': ColumnSpec('n_reads', 'read_count', ColumnType.INT, 0, None),
    'n_alns': ColumnSpec('n_alns', 'alignment_count', ColumnType.INT, 0, None),
    'n_reads_tad': ColumnSpec('n_reads_tad', 'read_count_tad', ColumnType.INT, 0, None),
    'n_bins': ColumnSpec('n_bins', 'bin_count', ColumnType.INT, 0, None),

    # Read Length Statistics
    'read_length_mean': ColumnSpec('read_length_mean', 'read_length_mean', ColumnType.FLOAT, 0, 100000),
    'read_length_std': ColumnSpec('read_length_std', 'read_length_std', ColumnType.FLOAT, 0, 100000),
    'read_length_min': ColumnSpec('read_length_min', 'read_length_min', ColumnType.INT, 0, 100000),
    'read_length_max': ColumnSpec('read_length_max', 'read_length_max', ColumnType.INT, 0, 100000),
    'read_length_median': ColumnSpec('read_length_median', 'read_length_median', ColumnType.FLOAT, 0, 100000),
    'read_length_mode': ColumnSpec('read_length_mode', 'read_length_mode', ColumnType.INT, 0, 100000),

    # GC Content (percentages)
    'gc_content_mean': ColumnSpec('gc_content_mean', 'gc_content_mean', ColumnType.FLOAT, 0, 100),
    'gc_content_std': ColumnSpec('gc_content_std', 'gc_content_std', ColumnType.FLOAT, 0, 100),
    'gc_content_total': ColumnSpec('gc_content_total', 'gc_content_total', ColumnType.FLOAT, 0, 100),

    # Quality Metrics
    'dust_mean': ColumnSpec('dust_mean', 'dust_mean', ColumnType.FLOAT, 0, None),
    'dust_std': ColumnSpec('dust_std', 'dust_std', ColumnType.FLOAT, 0, None),

    # Alignment Metrics
    'read_aligned_length': ColumnSpec('read_aligned_length', 'read_aligned_length_mean', ColumnType.FLOAT, 0, 100000),
    'read_aln_score': ColumnSpec('read_aln_score', 'read_alignment_score_mean', ColumnType.FLOAT, 0, None),
    'mapping_quality': ColumnSpec('mapping_quality', 'mapping_quality_mean', ColumnType.FLOAT, 0, 60),
    'edit_distances': ColumnSpec('edit_distances', 'edit_distance_mean', ColumnType.FLOAT, 0, None),

    # ANI (Average Nucleotide Identity) - percentages
    'read_ani_mean': ColumnSpec('read_ani_mean', 'read_ani_mean', ColumnType.FLOAT, 0, 100),
    'read_ani_std': ColumnSpec('read_ani_std', 'read_ani_std', ColumnType.FLOAT, 0, 100),
    'read_ani_median': ColumnSpec('read_ani_median', 'read_ani_median', ColumnType.FLOAT, 0, 100),

    # Coverage Metrics
    'bases_covered': ColumnSpec('bases_covered', 'bases_covered', ColumnType.INT, 0, None),
    'max_covered_bases': ColumnSpec('max_covered_bases', 'bases_covered_max', ColumnType.INT, 0, None),
    'mean_covered_bases': ColumnSpec('mean_covered_bases', 'bases_covered_mean', ColumnType.FLOAT, 0, None),
    'coverage_mean': ColumnSpec('coverage_mean', 'coverage_mean', ColumnType.FLOAT, 0, None),
    'coverage_mean_trunc': ColumnSpec('coverage_mean_trunc', 'coverage_mean_trimmed', ColumnType.FLOAT, 0, None),
    'coverage_mean_trunc_len': ColumnSpec('coverage_mean_trunc_len', 'coverage_mean_trimmed_length', ColumnType.INT, 0, None),
    'coverage_covered_mean': ColumnSpec('coverage_covered_mean', 'coverage_mean_covered_only', ColumnType.FLOAT, 0, None),

    # Reference Lengths
    'reference_length': ColumnSpec('reference_length', 'reference_length', ColumnType.INT, 0, None),
    'bam_reference_length': ColumnSpec('bam_reference_length', 'reference_length_bam', ColumnType.INT, 0, None),

    # Breadth Metrics (proportions 0-1)
    'breadth': ColumnSpec('breadth', 'breadth', ColumnType.FLOAT, 0, 1),
    'exp_breadth': ColumnSpec('exp_breadth', 'breadth_expected', ColumnType.FLOAT, 0, 1),
    'breadth_exp_ratio': ColumnSpec('breadth_exp_ratio', 'breadth_expected_ratio', ColumnType.FLOAT, 0, None),

    # Distribution Metrics
    'site_density': ColumnSpec('site_density', 'site_density', ColumnType.FLOAT, 0, None),
    'spatial_entropy': ColumnSpec('spatial_entropy', 'spatial_entropy', ColumnType.FLOAT, 0, None),
    'norm_spatial_entropy': ColumnSpec('norm_spatial_entropy', 'spatial_entropy_normalized', ColumnType.FLOAT, 0, 1),
    'gini': ColumnSpec('gini', 'gini_coefficient', ColumnType.FLOAT, 0, 1),
    'norm_gini': ColumnSpec('norm_gini', 'gini_coefficient_normalized', ColumnType.FLOAT, 0, 1),
    'c_v': ColumnSpec('c_v', 'coefficient_of_variation', ColumnType.FLOAT, 0, None),
    'd_i': ColumnSpec('d_i', 'diversity_index', ColumnType.FLOAT, 0, 1),
    'cov_evenness': ColumnSpec('cov_evenness', 'coverage_evenness', ColumnType.FLOAT, 0, 1),

    # Taxonomic Abundance
    'tax_abund_read': ColumnSpec('tax_abund_read', 'abundance_read_based', ColumnType.FLOAT, 0, 1),
    'tax_abund_aln': ColumnSpec('tax_abund_aln', 'abundance_alignment_based', ColumnType.FLOAT, 0, 1),
    'tax_abund_tad': ColumnSpec('tax_abund_tad', 'abundance_tad', ColumnType.FLOAT, 0, 1),
}

# Build reverse mapping: new_name -> old_name
NEW_TO_OLD_NAMES = {spec.new_name: old_name for old_name, spec in COLUMN_SPECS.items()}

# Get filterable columns
FILTERABLE_COLUMNS = {
    spec.new_name for spec in COLUMN_SPECS.values() if spec.filterable
}


class FilterSpec:
    """Represents a single column filter with min/max bounds."""

    def __init__(
        self,
        column: str,
        min_val: Optional[float],
        max_val: Optional[float],
        column_spec: ColumnSpec
    ):
        self.column = column
        self.min = min_val
        self.max = max_val
        self.column_spec = column_spec

    def __repr__(self):
        if self.min is not None and self.max is not None:
            return f"{self.column}: [{self.min}, {self.max}]"
        elif self.min is not None:
            return f"{self.column}: >= {self.min}"
        elif self.max is not None:
            return f"{self.column}: <= {self.max}"
        return f"{self.column}: no filter"


def get_column_index(column_name: str) -> int:
    """Get the column index in the RefStats struct for a given column name.

    Args:
        column_name: New column name

    Returns:
        Zero-based index of the column in RefStats struct (for get_column_value)

    Raises:
        KeyError: If column name is not found
    """
    # Map new column names to old internal names
    old_name = NEW_TO_OLD_NAMES.get(column_name)
    if old_name is None:
        raise KeyError(f"Unknown column: {column_name}")

    # This order MUST match the order in generic_filters.pyx get_column_value()
    STRUCT_COLUMN_ORDER = [
        'reference', 'n_reads', 'n_alns', 'read_length_mean', 'read_length_std',
        'read_length_min', 'read_length_max', 'read_length_median', 'read_length_mode',
        'gc_content_mean', 'gc_content_std', 'gc_content_total', 'dust_mean', 'dust_std',
        'read_aligned_length', 'read_aln_score', 'mapping_quality', 'edit_distances',
        'read_ani_mean', 'read_ani_std', 'read_ani_median',
        'bases_covered', 'max_covered_bases', 'mean_covered_bases',
        'coverage_mean', 'coverage_mean_trunc', 'coverage_mean_trunc_len', 'coverage_covered_mean',
        'reference_length', 'bam_reference_length',
        'breadth', 'exp_breadth', 'breadth_exp_ratio',
        'n_bins', 'site_density', 'spatial_entropy', 'norm_spatial_entropy', 'gini', 'norm_gini',
        'c_v', 'd_i', 'cov_evenness',
        'tax_abund_read', 'tax_abund_aln', 'tax_abund_tad', 'n_reads_tad'
    ]

    if old_name not in STRUCT_COLUMN_ORDER:
        raise KeyError(f"Column {old_name} not found in struct order")

    return STRUCT_COLUMN_ORDER.index(old_name)


def parse_filter_string(filter_str: str) -> Dict[str, FilterSpec]:
    """Parse filter specification string into structured filters.

    Args:
        filter_str: Comma-separated filters like "col1:10:100,col2:5:,col3::20"

    Returns:
        Dictionary mapping column name to FilterSpec

    Raises:
        ValueError: If filter syntax is invalid or column name unknown

    Examples:
        >>> parse_filter_string("read_ani_mean:90:")
        {'read_ani_mean': FilterSpec(read_ani_mean >= 90)}

        >>> parse_filter_string("breadth::0.95,coverage_mean:5:")
        {'breadth': FilterSpec(breadth <= 0.95), 'coverage_mean': FilterSpec(coverage_mean >= 5)}
    """
    filters = {}

    if not filter_str or not filter_str.strip():
        return filters

    for spec in filter_str.split(','):
        spec = spec.strip()
        if not spec:
            continue

        parts = spec.split(':')

        # Must have exactly 3 parts: column:min:max
        if len(parts) != 3:
            raise ValueError(
                f"Invalid filter syntax: '{spec}'\n"
                f"Expected format: 'column:min:max' (use empty for no bound)\n"
                f"Examples:\n"
                f"  'read_ani_mean:90:'      (>= 90)\n"
                f"  'breadth::0.95'          (<= 0.95)\n"
                f"  'coverage_mean:5:1000'   (between 5 and 1000)"
            )

        column, min_str, max_str = parts
        column = column.strip()

        # Validate column name exists and is filterable
        if column not in FILTERABLE_COLUMNS:
            from difflib import get_close_matches

            # Check if it's a known non-filterable column
            all_columns = {spec.new_name for spec in COLUMN_SPECS.values()}
            if column in all_columns:
                raise ValueError(
                    f"Column '{column}' is not filterable (it's a non-numeric identifier)"
                )

            # Suggest similar column names
            suggestions = get_close_matches(column, FILTERABLE_COLUMNS, n=3, cutoff=0.6)
            msg = f"Unknown column: '{column}'"
            if suggestions:
                msg += f"\nDid you mean: {', '.join(suggestions)}?"
            msg += f"\n\nAvailable filterable columns:\n"

            # Group by category for better readability
            categories = {
                'Read': ['read_count', 'read_length_mean', 'read_length_std', 'read_length_min',
                        'read_length_max', 'read_length_median', 'read_length_mode',
                        'read_ani_mean', 'read_ani_std', 'read_ani_median',
                        'read_aligned_length_mean', 'read_alignment_score_mean'],
                'Coverage': ['coverage_mean', 'coverage_mean_trimmed', 'coverage_mean_covered_only',
                           'coverage_evenness', 'breadth', 'breadth_expected', 'breadth_expected_ratio',
                           'bases_covered', 'bases_covered_max', 'bases_covered_mean'],
                'Quality': ['mapping_quality_mean', 'edit_distance_mean', 'gc_content_mean',
                          'gc_content_std', 'dust_mean', 'dust_std'],
                'Distribution': ['spatial_entropy', 'spatial_entropy_normalized', 'gini_coefficient',
                               'gini_coefficient_normalized', 'coefficient_of_variation',
                               'diversity_index', 'site_density'],
                'Other': ['alignment_count', 'reference_length', 'reference_length_bam',
                         'bin_count', 'abundance_read_based', 'abundance_alignment_based',
                         'abundance_tad', 'read_count_tad']
            }

            for cat, cols in categories.items():
                matching = [c for c in cols if c in FILTERABLE_COLUMNS]
                if matching:
                    msg += f"\n  {cat}: {', '.join(matching)}"

            raise ValueError(msg)

        # Get column spec
        old_name = NEW_TO_OLD_NAMES[column]
        col_spec = COLUMN_SPECS[old_name]

        # Parse bounds
        min_val = None
        max_val = None

        min_str = min_str.strip()
        max_str = max_str.strip()

        if min_str:
            try:
                min_val = float(min_str) if col_spec.dtype == ColumnType.FLOAT else int(min_str)
            except ValueError:
                raise ValueError(
                    f"Invalid minimum value '{min_str}' for column '{column}'\n"
                    f"Expected {col_spec.dtype.value} type"
                )

        if max_str:
            try:
                max_val = float(max_str) if col_spec.dtype == ColumnType.FLOAT else int(max_str)
            except ValueError:
                raise ValueError(
                    f"Invalid maximum value '{max_str}' for column '{column}'\n"
                    f"Expected {col_spec.dtype.value} type"
                )

        # At least one bound must be specified
        if min_val is None and max_val is None:
            raise ValueError(
                f"Filter for '{column}' has no bounds specified\n"
                f"Use 'column:min:' or 'column::max' or 'column:min:max'"
            )

        # Validate min < max
        if min_val is not None and max_val is not None:
            if min_val >= max_val:
                raise ValueError(
                    f"Invalid range for '{column}': min ({min_val}) must be < max ({max_val})"
                )

        # Validate against column constraints
        if min_val is not None:
            if col_spec.min_value is not None and min_val < col_spec.min_value:
                raise ValueError(
                    f"Minimum value {min_val} for '{column}' is below valid range "
                    f"(must be >= {col_spec.min_value})"
                )
            if col_spec.max_value is not None and min_val > col_spec.max_value:
                raise ValueError(
                    f"Minimum value {min_val} for '{column}' is above valid range "
                    f"(must be <= {col_spec.max_value})"
                )

        if max_val is not None:
            if col_spec.min_value is not None and max_val < col_spec.min_value:
                raise ValueError(
                    f"Maximum value {max_val} for '{column}' is below valid range "
                    f"(must be >= {col_spec.min_value})"
                )
            if col_spec.max_value is not None and max_val > col_spec.max_value:
                raise ValueError(
                    f"Maximum value {max_val} for '{column}' is above valid range "
                    f"(must be <= {col_spec.max_value})"
                )

        filters[column] = FilterSpec(column, min_val, max_val, col_spec)

    return filters


def get_tsv_header() -> str:
    """Generate TSV header line with new column names.

    Returns:
        Tab-separated header string
    """
    return '\t'.join(spec.new_name for spec in COLUMN_SPECS.values())


def export_filters_to_cython_format(filters: Dict[str, FilterSpec]) -> list:
    """Convert filter dict to format suitable for Cython processing.

    Returns a list of tuples: [(column_index, min_value, max_value), ...]
    where column_index is the position in the TSV output.

    Args:
        filters: Dictionary of FilterSpec objects keyed by column name

    Returns:
        List of (column_index, min_value, max_value) tuples
    """
    result = []
    for column_name, filter_spec in filters.items():
        col_idx = get_column_index(column_name)
        min_val = filter_spec.min if filter_spec.min is not None else float('-inf')
        max_val = filter_spec.max if filter_spec.max is not None else float('inf')
        result.append((col_idx, min_val, max_val))

    return result
