"""Read reassignment module for filterBAM.

This module provides the command-line interface for reassigning multi-mapping reads
using an Expectation-Maximization (EM) algorithm with optional SQUAREM acceleration.
"""

import os
import sys
import logging
from time import perf_counter
from typing import Any, Dict, List, Optional, Tuple

from bam_filter.utils import (
    get_arguments,
    create_output_files,
    handle_warning,
    is_debug,
)
from bam_filter.processor import process_bam_with_em
from bam_filter import logging as bf_logging

LOG_TAG = "REASSIGN"


def _info(message: str, *args: Any) -> None:
    bf_logging.log(LOG_TAG, message, *args)


def _warn(message: str, *args: Any) -> None:
    bf_logging.warn(message, *args)


def _error(message: str, *args: Any) -> None:
    bf_logging.error(message, *args)


def _debug(level: int, message: str, *args: Any) -> None:
    bf_logging.verbose(level, LOG_TAG, message, *args)


log = logging.getLogger("bam_filter")


def reassign_reads(
    bam_file: str,
    output_bam: Optional[str] = None,
    num_threads: int = 1,
    verbose: bool = False,
    # PMD parameters
    calculate_pmd: bool = True,
    library_type: str = "ds",
    # Read filtering parameters
    min_read_count: int = 1,
    min_read_ani: float = 0.0,
    min_read_length: int = 30,
    max_read_length: int = 10000,
    # EM algorithm parameters
    max_em_iterations: int = 25,
    em_tolerance: float = 1e-5,
    min_probability: float = 1e-6,
    prob_fraction: float = 0.1,
    prior_weight: float = 0.01,
    use_squarem_acceleration: bool = True,
    # SQUAREM parameters
    enable_globalization: bool = True,
    squarem_start_iter: int = 2,
    backtrack_factor: float = 0.5,
    max_backtrack_steps: int = 5,
    steplength_scheme: int = 3,
    # Reference length override
    reference_lengths_tsv: Optional[str] = None,
    reference_stats_tsv: Optional[str] = None,
    information_threshold: float = -999.0,
    # Graph construction parameters
    graph_min_edge_weight: int = 0,
    # Clustering parameters
    clustering: bool = False,
    community_resolution: float = 1.0,
    community_max_iterations: int = 10,
    outlier_method: str = "mad",
    # Graph export
    graph_export: Optional[str] = None,
    # Taxonomy parameters
    taxonomy_db: Optional[str] = None,
    taxonomy_min_rank: int = 6,
    taxonomy_cross_domain_threshold: float = 0.10,
    taxonomy_kingdom_threshold: float = 0.25,
    taxonomy_genus_threshold: float = 0.50,
    # Taxonomy-informed filtering parameters
    taxonomy_filter: bool = False,
    taxonomy_strict_filter: bool = True,
    taxonomy_strict_min_connections: int = 5,
    taxonomy_weighted_outlier: bool = True,
    taxonomy_anomaly_weight: float = 2.0,
    taxonomy_second_chance: bool = True,
    taxonomy_second_chance_cc: float = 0.3,
    # Cross-domain removal options (None = auto, True = enabled, False = disabled)
    remove_cross_domain_alignments: bool = None,
    remove_cross_domain_references: bool = None,
    remove_cross_domain_all: bool = None,
    no_cross_domain_removal: bool = False,
    detect_misannotations: bool = None,
) -> Dict[str, Any]:
    # Handle --remove-cross-domain-all as shorthand for both
    if remove_cross_domain_all:
        remove_cross_domain_alignments = True
        remove_cross_domain_references = True

    # Auto-enable combined mode when taxonomy filtering is active with a database
    # Unless explicitly disabled with --no-cross-domain-removal
    if taxonomy_filter and taxonomy_db and not no_cross_domain_removal:
        # If no explicit cross-domain flags set, enable combined mode by default
        if remove_cross_domain_alignments is None and remove_cross_domain_references is None:
            remove_cross_domain_alignments = True
            remove_cross_domain_references = True
            detect_misannotations = True if detect_misannotations is None else detect_misannotations
            _info(
                "Cross-domain removal automatically enabled (combined mode). "
                "Use --no-cross-domain-removal to disable."
            )

    # Convert None to False for downstream code
    if remove_cross_domain_alignments is None:
        remove_cross_domain_alignments = False
    if remove_cross_domain_references is None:
        remove_cross_domain_references = False
    if detect_misannotations is None:
        detect_misannotations = False

    # Validate dependencies
    if (remove_cross_domain_alignments or remove_cross_domain_references or detect_misannotations):
        if not taxonomy_db:
            raise ValueError(
                "Cross-domain removal flags require a taxonomy database (--taxonomy-db)."
            )
        if not taxonomy_filter:
            taxonomy_filter = True
            _info(
                "Taxonomy filtering automatically enabled because "
                "cross-domain removal was requested."
            )

    """Reassign multi-mapping reads using EM algorithm with SQUAREM acceleration.

    This function wraps the high-performance C/Cython processor that handles
    BAM file reading, alignment scoring, EM algorithm execution, and result writing.

    Parameters
    ----------
    bam_file : str
        Path to input BAM file
    output_bam : str, optional
        Path to output BAM file
    num_threads : int, default=1
        Number of parallel threads
    verbose : bool, default=False
        Enable detailed logging
    calculate_pmd : bool, default=True
        Calculate Post-Mortem Damage scores for ancient DNA
    library_type : str, default="ds"
        Library type: "ds" (double-stranded) or "ss" (single-stranded)
    min_read_count : int, default=1
        Minimum reads per reference
    min_read_ani : float, default=0.0
        Minimum average nucleotide identity (%)
    min_read_length : int, default=30
        Minimum read length (bp)
    max_read_length : int, default=10000
        Maximum read length (bp)
    max_em_iterations : int, default=25
        Maximum EM iterations
    em_tolerance : float, default=1e-5
        Convergence tolerance
    min_probability : float, default=1e-6
        Minimum probability threshold
    prob_fraction : float, default=0.1
        Fraction of alignments to keep by probability
    prior_weight : float, default=0.01
        Prior weight for EM initialization
    use_squarem_acceleration : bool, default=True
        Enable SQUAREM acceleration (Varadhan & Roland 2008)
    enable_globalization : bool, default=True
        Enable globalization with backtracking
    squarem_start_iter : int, default=2
        Iteration to start SQUAREM
    backtrack_factor : float, default=0.5
        Backtracking step size factor
    max_backtrack_steps : int, default=5
        Maximum backtracking iterations
    steplength_scheme : int, default=3
        Steplength scheme (1=S1, 2=S2, 3=S3 recommended)
    reference_lengths_tsv : str, optional
        TSV file with reference lengths (reference_name<tab>length)
    reference_stats_tsv : str, optional
        Output path for reference statistics TSV
    graph_min_edge_weight : int, default=0
        Minimum edge weight for graph (0=auto via elbow detection, -1=none, >0=explicit)
    clustering : bool, default=False
        Enable Community clustering
    community_resolution : float, default=1.0
        Community resolution parameter
    community_max_iterations : int, default=10
        Maximum Community iterations
    graph_export : str, optional
        Path to export graph in GraphML format

    Returns
    -------
    dict
        Processing results containing alignment counts, EM statistics,
        and processing metadata
    """
    overall_start = perf_counter()

    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"Input BAM file not found: {bam_file}")

    if reference_lengths_tsv and not os.path.exists(reference_lengths_tsv):
        _warn("Reference lengths TSV not found: %s", reference_lengths_tsv)
        reference_lengths_tsv = None

    output_display = output_bam or "not specified"
    _info(
        "Processing configuration: %s → %s",
        bam_file,
        output_display,
    )
    _info(
        "Read quality filters: minimum coverage=%d, minimum identity=%.1f%%, length range=%d-%d bp",
        min_read_count,
        min_read_ani,
        min_read_length,
        max_read_length,
    )
    _info(
        "EM algorithm settings: probability threshold=%.1e, retention fraction=%.2f, prior weight=%.2f, SQUAREM acceleration=%s, globalization=%s, step scheme=S%d",
        min_probability,
        prob_fraction,
        prior_weight,
        "enabled" if use_squarem_acceleration else "disabled",
        "enabled" if enable_globalization else "disabled",
        steplength_scheme,
    )

    try:
        processing_start = perf_counter()
        pipeline_steps = [
            "alignment scoring",
            "memory optimization",
            "EM optimization",
            "probability filtering",
        ]
        if reference_stats_tsv:
            pipeline_steps.append("connectivity analysis")
        if clustering:
            pipeline_steps.append("cluster-aware refinement")
        if output_bam:
            pipeline_steps.append("BAM file generation")

        _info("Pipeline: %s", " → ".join(pipeline_steps))

        taxdb_path = None
        accmap_path = None
        if taxonomy_db:
            if not os.path.exists(taxonomy_db):
                _warn("Taxonomy database path does not exist: %s", taxonomy_db)
                _warn("Continuing without taxonomy-aware graph analysis")
            else:
                accession_map_path = os.path.join(taxonomy_db, "accession_map.parquet")
                if os.path.exists(accession_map_path):
                    taxdb_path = taxonomy_db
                    accmap_path = accession_map_path
                else:
                    _warn("No accession_map.parquet found in %s", taxonomy_db)
                    _warn("Continuing without taxonomy-aware graph analysis")

        result = process_bam_with_em(
            bam_file=bam_file,
            output_bam=output_bam,
            num_threads=num_threads,
            verbose=verbose,
            calculate_pmd=calculate_pmd,
            library_type=library_type,
            min_read_count=min_read_count,
            min_read_length=min_read_length,
            max_read_length=max_read_length,
            min_read_ani=min_read_ani,
            max_em_iterations=max_em_iterations,
            em_tolerance=em_tolerance,
            min_probability=min_probability,
            prob_fraction=prob_fraction,
            prior_weight=prior_weight,
            use_squarem_acceleration=use_squarem_acceleration,
            enable_globalization=enable_globalization,
            squarem_start_iter=squarem_start_iter,
            backtrack_factor=backtrack_factor,
            max_backtrack_steps=max_backtrack_steps,
            steplength_scheme=steplength_scheme,
            reference_lengths_tsv=reference_lengths_tsv,
            reference_stats_tsv=reference_stats_tsv,
            init_prior_strength=0.1,
            information_threshold=information_threshold,
            graph_min_edge_weight=graph_min_edge_weight,
            clustering=clustering,
            community_resolution=community_resolution,
            community_max_iterations=community_max_iterations,
            outlier_method=outlier_method,
            graph_export=graph_export,
            taxonomy_db=taxdb_path,
            taxonomy_accession_map=accmap_path,
            taxonomy_min_rank=taxonomy_min_rank,
            taxonomy_cross_domain_threshold=taxonomy_cross_domain_threshold,
            taxonomy_kingdom_threshold=taxonomy_kingdom_threshold,
            taxonomy_genus_threshold=taxonomy_genus_threshold,
            taxonomy_filter_enabled=taxonomy_filter,
            taxonomy_strict_filter=taxonomy_strict_filter,
            taxonomy_strict_min_connections=taxonomy_strict_min_connections,
            taxonomy_weighted_outlier=taxonomy_weighted_outlier,
            taxonomy_anomaly_weight=taxonomy_anomaly_weight,
            taxonomy_second_chance=taxonomy_second_chance,
            taxonomy_second_chance_cc=taxonomy_second_chance_cc,
            remove_cross_domain_alignments=remove_cross_domain_alignments,
            remove_cross_domain_references=remove_cross_domain_references,
            detect_misannotations=detect_misannotations,
        )
        processing_duration = perf_counter() - processing_start
        _info("Core processing completed in %.2f seconds", processing_duration)

        if not result.get("success", False):
            error_msg = result.get("error_message")
            if error_msg is None:
                error_msg_str = "Unknown error in unified processor"
            else:
                error_msg_str = str(error_msg)

            _error("Analysis pipeline failed: %s", error_msg_str)
            _debug(1, "Full result: %r", result)
            raise RuntimeError(f"BAM processing failed: {error_msg_str}")

        overall_duration = perf_counter() - overall_start
        _info("Total execution time: %.2f seconds", overall_duration)

        return {
            "success": result.get("success", False),
            "n_processed_alignments": result.get("n_input_alignments", 0),
            "n_filtered_alignments": result.get("n_filtered_alignments", 0),
            "n_unique_reads": result.get("n_unique_reads", 0),
            "n_final_unique_reads": result.get("n_final_unique_reads", 0),
            "n_unique_references": result.get("n_unique_references", 0),
            "em_iterations": result.get("em_iterations", 0),
            "em_converged": result.get("em_converged", False),
            "final_likelihood": result.get("final_likelihood", 0.0),
            "processing_time": result.get("total_time", 0.0),
            "peak_memory_mb": 0,
            "squarem_used": result.get("squarem_used", False),
            "globalization_enabled": result.get("globalization_enabled", False),
            "steplength_scheme": result.get("steplength_scheme", "Unknown"),
            "paper_aligned": result.get("paper_aligned", False),
            "implementation": result.get("implementation", "Unknown"),
            "pmd_enabled": result.get("pmd_enabled", False),
            "library_type": result.get("library_type", "ds"),
            "pmd_storage_mb": result.get("pmd_storage_mb", 0.0),
            "pmd_tags_written": result.get("pmd_tags_written", False),
            "reference_lengths_tsv": reference_lengths_tsv,
            "processing_duration": processing_duration,
            "overall_duration": overall_duration,
        }

    except Exception as e:
        _error("BAM processing failed: %s", e)
        raise


def parse_memory_string(mem_str: str) -> int:
    """Parse memory string (e.g. '2G', '500M') into bytes.

    Parameters
    ----------
    mem_str : str
        Memory string with optional suffix (G, M, K)

    Returns
    -------
    int
        Memory size in bytes
    """
    if mem_str is None:
        return None

    mem_str = str(mem_str).strip().upper()
    if mem_str.endswith("G"):
        return int(float(mem_str[:-1]) * 1024**3)
    elif mem_str.endswith("M"):
        return int(float(mem_str[:-1]) * 1024**2)
    elif mem_str.endswith("K"):
        return int(float(mem_str[:-1]) * 1024)
    else:
        return int(mem_str)


def reassign(args):
    """Main entry point for read reassignment command.

    This function provides the command-line interface, parsing arguments
    and delegating processing to reassign_reads().

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments
    """
    _info("Starting read reassignment analysis")

    # Parse input file path
    if hasattr(args, "bam"):
        bam_file = args.bam
    elif hasattr(args, "bam_file"):
        bam_file = args.bam_file
    else:
        bam_file = args.input

    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"Input BAM file not found: {bam_file}")

    # Determine output BAM file
    output_bam = None
    if hasattr(args, "bam_reassigned") and args.bam_reassigned:
        output_bam = args.bam_reassigned
    elif hasattr(args, "output_bam") and args.output_bam:
        output_bam = args.output_bam
    else:
        base_name = os.path.splitext(bam_file)[0]
        output_bam = f"{base_name}_reassigned.bam"

    taxonomy_filter_enabled = getattr(args, "taxonomy_filter_enabled", False)
    taxonomy_strict_filter = getattr(args, "taxonomy_strict_filter", None)
    if taxonomy_strict_filter is None:
        taxonomy_strict_filter = taxonomy_filter_enabled

    taxonomy_weighted_outlier = getattr(args, "taxonomy_weighted_outlier", None)
    if taxonomy_weighted_outlier is None:
        taxonomy_weighted_outlier = taxonomy_filter_enabled

    taxonomy_second_chance = getattr(args, "taxonomy_second_chance", None)
    if taxonomy_second_chance is None:
        taxonomy_second_chance = taxonomy_filter_enabled

    # Extract parameters from args
    params = {
        "bam_file": bam_file,
        "output_bam": output_bam,
        "num_threads": getattr(args, "num_threads", getattr(args, "threads", 1)),
        "verbose": is_debug(),
        "calculate_pmd": not getattr(args, "disable_pmd", False),
        "library_type": getattr(args, "library_type", "ds"),
        "min_read_count": getattr(args, "min_read_count", 1),
        "min_read_ani": getattr(args, "min_read_ani", 0.0),
        "min_read_length": getattr(args, "min_read_length", 30),
        "max_read_length": getattr(args, "max_read_length", 10000),
        "max_em_iterations": getattr(
            args, "max_em_iterations", getattr(args, "reassign_iters", 25)
        ),
        "em_tolerance": getattr(args, "em_tolerance", 1e-5),
        "min_probability": getattr(
            args, "min_probability", getattr(args, "min_prob", 1e-6)
        ),
        "prob_fraction": getattr(args, "prob_fraction", 0.1),
        "prior_weight": getattr(args, "prior_weight", 0.01),
        "use_squarem_acceleration": getattr(
            args,
            "use_squarem_acceleration",
            getattr(args, "acceleration", "squarem") == "squarem",
        ),
        "enable_globalization": getattr(args, "enable_globalization", True),
        "squarem_start_iter": getattr(args, "squarem_start_iter", 2),
        "backtrack_factor": getattr(args, "backtrack_factor", 0.5),
        "max_backtrack_steps": getattr(args, "max_backtrack_steps", 5),
        "steplength_scheme": getattr(args, "steplength_scheme", 3),
        "reference_lengths_tsv": getattr(args, "reference_lengths_tsv", None),
        "reference_stats_tsv": getattr(args, "reference_stats_tsv", None),
        "graph_min_edge_weight": getattr(args, "graph_min_edge_weight", 0),
        "clustering": getattr(args, "clustering", False),
        "community_resolution": getattr(args, "community_resolution", 1.0),
        "community_max_iterations": getattr(args, "community_max_iterations", 10),
        "outlier_method": getattr(args, "outlier_method", "mad"),
        "graph_export": getattr(args, "graph_export", None),
        # Taxonomy parameters
        "taxonomy_db": getattr(args, "taxonomy_db", None),
        "taxonomy_min_rank": getattr(args, "taxonomy_min_rank", 6),
        "taxonomy_cross_domain_threshold": getattr(args, "taxonomy_cross_domain_threshold", 0.10),
        "taxonomy_kingdom_threshold": getattr(args, "taxonomy_kingdom_threshold", 0.25),
        "taxonomy_genus_threshold": getattr(args, "taxonomy_genus_threshold", 0.50),
        # Taxonomy-informed filtering parameters
        "taxonomy_filter": taxonomy_filter_enabled,
        "taxonomy_strict_filter": taxonomy_strict_filter,
        "taxonomy_strict_min_connections": getattr(args, "taxonomy_strict_min_connections", 5),
        "taxonomy_weighted_outlier": taxonomy_weighted_outlier,
        "taxonomy_anomaly_weight": getattr(args, "taxonomy_anomaly_weight", 2.0),
        "taxonomy_second_chance": taxonomy_second_chance,
        "taxonomy_second_chance_cc": getattr(args, "taxonomy_second_chance_cc", 0.3),
        # Cross-domain removal options
        "remove_cross_domain_alignments": getattr(args, "remove_cross_domain_alignments", False),
        "remove_cross_domain_references": getattr(args, "remove_cross_domain_references", False),
        "remove_cross_domain_all": getattr(args, "remove_cross_domain_all", False),
        "detect_misannotations": getattr(args, "detect_misannotations", False),
    }

    try:
        result = reassign_reads(**params)

        if result.get("success"):
            _info("Read reassignment completed successfully")
        else:
            _warn("Read reassignment finished with warnings")

        _info("Output file: %s", output_bam)
        return result

    except Exception as e:
        _error("Read reassignment analysis failed: %s", e)
        raise
