"""
Helper utilities for running the Cython LCA stats processor from Python code.
"""

from bam_filter import logging as bf_logging
from bam_filter.processor_lca_stats import (
    process_lca_stats_wrapper,
    configure_lca_stats_thresholds,
)


def run_lca_stats(
    bam_path: str,
    lca_per_read_path: str,
    output_path: str,
    taxonomy_db_path: str,
    num_threads: int = 1,
    verbose: bool = False,
    min_read_ani: float = 0.0,
    min_read_length: int = 0,
    max_read_length: int = (2**31) - 1,
    scale: int = 1_000_000,
    trim_ends: bool = False,
    trim_min: int = 10,
    trim_max: int = 90,
    taxdb = None,
):
    from bam_filter.taxonomy_db import TaxonomyDatabase

    if taxdb is None:
        if verbose:
            bf_logging.log("LCA_STATS", f"Loading taxonomy from {taxonomy_db_path}")

        # Load taxonomy database from parquet
        taxdb = TaxonomyDatabase.from_parquet(taxonomy_db_path)
    else:
        if verbose:
            bf_logging.log("LCA_STATS", "Reusing pre-loaded taxonomy database")

    if verbose:
        bf_logging.log("LCA_STATS", f"Processing BAM: {bam_path}")
        bf_logging.log("LCA_STATS", f"LCA assignments: {lca_per_read_path}")
        bf_logging.log("LCA_STATS", f"Output: {output_path}")

    configure_lca_stats_thresholds(
        min_read_ani=min_read_ani,
        min_read_length=min_read_length,
        max_read_length=max_read_length,
        scale=scale,
        trim_ends=int(trim_ends),
        trim_min=trim_min,
        trim_max=trim_max,
    )

    # Run processing using the proven C implementation
    process_lca_stats_wrapper(
        bam_path.encode("utf-8"),
        lca_per_read_path.encode("utf-8"),
        output_path.encode("utf-8"),
        taxdb,
        taxonomy_db_path.encode("utf-8"),
        num_threads,
        verbose,
    )
