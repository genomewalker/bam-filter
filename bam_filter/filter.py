"""Reference filtering module for filterBAM.

This module provides the command-line interface for calculating coverage statistics
and filtering references based on evenness criteria.
"""

from __future__ import annotations

import os
import sys
from typing import Any

from bam_filter import logging as bf_logging
from bam_filter import stats
from bam_filter.utils import defaults
from bam_filter.filter_parser import (
    parse_filter_string,
    FILTERABLE_COLUMNS,
    COLUMN_SPECS,
    export_filters_to_cython_format,
)

LOG_TAG = "FILTER"


def _info(message: str, *args: Any) -> None:
    bf_logging.log(LOG_TAG, message, *args)


def _warn(message: str, *args: Any) -> None:
    bf_logging.warn(message, *args)


def _error(message: str, *args: Any) -> None:
    bf_logging.error(message, *args)


def _debug(level: int, message: str, *args: Any) -> None:
    bf_logging.verbose(level, LOG_TAG, message, *args)


def list_filterable_columns():
    """Print all available filterable columns and their specifications."""
    print("\n=== Filterable Columns ===\n")
    print("Filter Name → TSV Column Name [type] (constraints)\n")

    # Group by category
    categories = {
        'Read Statistics': [
            'read_count', 'alignment_count', 'read_length_mean', 'read_length_std',
            'read_length_min', 'read_length_max', 'read_length_median', 'read_length_mode',
            'read_ani_mean', 'read_ani_std', 'read_ani_median',
            'read_aligned_length_mean', 'read_alignment_score_mean'
        ],
        'Coverage Metrics': [
            'coverage_mean', 'coverage_mean_trimmed', 'coverage_mean_covered_only',
            'coverage_evenness', 'breadth', 'breadth_expected', 'breadth_expected_ratio',
            'bases_covered', 'bases_covered_max', 'bases_covered_mean'
        ],
        'Quality Metrics': [
            'mapping_quality_mean', 'edit_distance_mean', 'gc_content_mean',
            'gc_content_std', 'gc_content_total', 'dust_mean', 'dust_std'
        ],
        'Distribution Metrics': [
            'spatial_entropy', 'spatial_entropy_normalized', 'gini_coefficient',
            'gini_coefficient_normalized', 'coefficient_of_variation',
            'diversity_index', 'site_density'
        ],
        'Reference & Other': [
            'reference_length', 'reference_length_bam', 'bin_count',
            'abundance_read_based', 'abundance_alignment_based',
            'abundance_tad', 'read_count_tad'
        ]
    }

    for category, columns in categories.items():
        print(f"\n{category}:")
        for col in columns:
            if col in FILTERABLE_COLUMNS:
                # Find the column spec
                for old_name, spec in COLUMN_SPECS.items():
                    if spec.new_name == col:
                        bounds = []
                        if spec.min_value is not None:
                            bounds.append(f"min={spec.min_value}")
                        if spec.max_value is not None:
                            bounds.append(f"max={spec.max_value}")
                        bounds_str = f" ({', '.join(bounds)})" if bounds else ""
                        # Show: filter_name → tsv_column_name [type] (constraints)
                        print(f"  {col:35} → {old_name:30} [{spec.dtype.value}]{bounds_str}")
                        break

    print("\n\nNotes:")
    print("  - Use the filter name (left side) in --filter arguments")
    print("  - TSV output files use the column name (right side)")
    print("\nExample usage:")
    print("  --filter 'read_ani_mean:90:'")
    print("  --filter 'breadth::0.95,coverage_mean:5:'")
    print("  --filter 'read_count:10:,coverage_evenness:0.1:0.9'\n")


def filter_references(args):
    """Main entrypoint for the filter subcommand.

    Calculates comprehensive per-reference statistics and applies filtering
    criteria to identify references with acceptable coverage patterns.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments from get_arguments()

    Returns
    -------
    dict or int
        Processing results containing statistics and filter outcomes

    Raises
    ------
    FileNotFoundError
        If input BAM file does not exist
    RuntimeError
        If filtering fails
    """
    # Handle --list-columns (doesn't need BAM file)
    if getattr(args, "list_columns", False):
        list_filterable_columns()
        sys.exit(0)

    # Determine BAM file path
    bam_file = getattr(args, "bam", None)
    if bam_file is None and hasattr(args, "bam_file"):
        bam_file = getattr(args, "bam_file")

    # Validate input exists (required for all operations except --list-columns)
    if bam_file is None:
        raise ValueError("--bam argument is required (unless using --list-columns)")
    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"Input BAM file not found: {bam_file}")

    verbosity = getattr(args, "verbosity", bf_logging.get_level())
    bf_logging.set_level(verbosity)
    verbose_flag = verbosity >= bf_logging.LogLevel.INFO

    output_targets = [
        path
        for path in (
            getattr(args, "output", None),
            getattr(args, "filtered_output", None),
            getattr(args, "filtered_bam", None),
        )
        if path
    ]

    _info("Starting reference filtering for %s", bam_file)
    _info(
        "Outputs: %s",
        ", ".join(output_targets) if output_targets else "none (stats only)",
    )
    _debug(1, "Verbosity level set to %s", bf_logging.LogLevel(verbosity).name)
    pipeline_steps = ["input validation", "statistics aggregation", "filter evaluation"]
    if getattr(args, "filtered_bam", None):
        pipeline_steps.append("filtered BAM writing")
    _info("Pipeline: %s", " → ".join(pipeline_steps))

    # Parse generic filter specification
    filter_spec_str = getattr(args, "filter_spec", None)
    min_read_count = getattr(args, "min_read_count", None)
    generic_filters_list = None
    filters_dict = {}

    # Handle --filter specification
    if filter_spec_str:
        try:
            _info("Parsing filter specification: %s", filter_spec_str)
            filters_dict = parse_filter_string(filter_spec_str)
        except ValueError as e:
            _error("Filter parsing error: %s", str(e))
            raise RuntimeError(f"Invalid filter specification: {e}") from e

    # Note: min_read_count is passed to compute_bam_stats for index-based pre-filtering
    # It is NOT added to generic_filters because it's applied earlier (before stats computation)
    # to skip references with too few reads using the BAM index
    if min_read_count is not None and min_read_count > 0:
        _info("Using -n/--min-read-count for index-based pre-filtering: %d", min_read_count)

    # Convert to Cython format
    if filters_dict:
        generic_filters_list = export_filters_to_cython_format(filters_dict)
        _info("Active filter(s):")
        for col_name, filt in filters_dict.items():
            _info("  - %s", str(filt))

    # Call Cython stats module
    result = stats.compute_bam_stats(
        bam_file=bam_file,
        batch_size_param=getattr(args, "batch_size", 100),
        verbose=verbose_flag,
        verbosity_level=int(verbosity),
        num_threads=getattr(args, "threads", getattr(args, "num_threads", 1)),
        show_progress=getattr(args, "show_progress", False),
        min_read_length=getattr(args, "min_read_length", defaults["min_read_length"]),
        max_read_length=getattr(args, "max_read_length", defaults["max_read_length"]),
        min_read_ani=getattr(args, "min_read_ani", defaults["min_read_ani"]),
        min_read_count=min_read_count if min_read_count is not None else 1,
        generic_filters=generic_filters_list,  # New parameter
        output=getattr(args, "output", None),
        filtered_output=getattr(args, "filtered_output", None),
        filtered_bam=getattr(args, "filtered_bam", None),
        scale=getattr(args, "scale", defaults["scale"]),
        trim_ends=getattr(args, "trim_ends", 0),
        trim_min=getattr(args, "trim_min", 10),
        trim_max=getattr(args, "trim_max", 90),
        reference_lengths_tsv=getattr(args, "reference_lengths_tsv", None),
    )

    # Check result and raise on failure
    if isinstance(result, dict):
        if not result.get("success", True):
            raise RuntimeError(
                f"Filtering failed: {result.get('error_message', 'unknown')}"
            )
        _info(
            "Filtering completed: %d references processed",
            result.get("n_total_references", 0),
        )
        return result

    if isinstance(result, int) and result != 0:
        raise RuntimeError(f"Filtering returned non-zero status: {result}")

    _info("Filtering completed")
    return result
