"""Reference filtering module for filterBAM.

This module provides the command-line interface for calculating coverage statistics
and filtering references based on evenness criteria.
"""

from __future__ import annotations

import os
from typing import Any

from bam_filter import logging as bf_logging
from bam_filter import stats
from bam_filter.utils import defaults

LOG_TAG = "FILTER"


def _info(message: str, *args: Any) -> None:
    bf_logging.log(LOG_TAG, message, *args)


def _warn(message: str, *args: Any) -> None:
    bf_logging.warn(message, *args)


def _error(message: str, *args: Any) -> None:
    bf_logging.error(message, *args)


def _debug(level: int, message: str, *args: Any) -> None:
    bf_logging.verbose(level, LOG_TAG, message, *args)


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
    # Determine BAM file path
    bam_file = getattr(args, "bam", None)
    if bam_file is None and hasattr(args, "bam_file"):
        bam_file = getattr(args, "bam_file")

    # Validate input exists
    if bam_file is None or not os.path.exists(bam_file):
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
        min_read_count=getattr(args, "min_read_count", defaults["min_read_count"]),
        filter_min_avg_read_ani=getattr(args, "min_avg_read_ani", None),
        filter_min_expected_breadth_ratio=getattr(
            args, "min_expected_breadth_ratio", None
        ),
        filter_min_breadth=getattr(args, "min_breadth", None),
        filter_min_coverage_evenness=getattr(args, "min_coverage_evenness", None),
        filter_max_coeff_var=getattr(args, "min_coeff_var", None),
        filter_min_coverage_mean=getattr(args, "min_coverage_mean", None),
        filter_min_norm_entropy=getattr(args, "min_norm_entropy", None),
        filter_max_norm_gini=getattr(args, "min_norm_gini", None),
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
        return result

    # Handle non-dict return (legacy compatibility)
    if isinstance(result, int) and result != 0:
        raise RuntimeError(f"Filtering returned non-zero status: {result}")

    if isinstance(result, dict):
        _info(
            "Filtering completed: %d references processed",
            result.get("n_total_references", 0),
        )
    else:
        _info("Filtering completed")

    return result
