import argparse
import sys
import gzip
import os
import shutil
# import pandas as pd
from multiprocessing import Pool
from bam_filter import logging as bf_logging
from functools import partial
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from os import devnull
import tqdm
from bam_filter import __version__
import time
from itertools import chain
import numpy as np
from pathlib import Path
import tempfile
from difflib import get_close_matches
from typing import Optional

LOG_TAG = "UTILS"

_VERBOSITY_NAME_TO_LEVEL = {
    "quiet": bf_logging.LogLevel.QUIET,
    "summary": bf_logging.LogLevel.SUMMARY,
    "info": bf_logging.LogLevel.INFO,
    "debug": bf_logging.LogLevel.DEBUG,
    "trace": bf_logging.LogLevel.TRACE,
}

_VERBOSITY_COUNT_TO_LEVEL = {
    0: bf_logging.LogLevel.SUMMARY,
    1: bf_logging.LogLevel.INFO,
    2: bf_logging.LogLevel.DEBUG,
}


def _info(message: str) -> None:
    bf_logging.log(LOG_TAG, message)


def _warn(message: str) -> None:
    bf_logging.warn(f"{LOG_TAG}: {message}")


def _error(message: str) -> None:
    bf_logging.error(f"{LOG_TAG}: {message}")


def handle_warning(message, category, filename, lineno, file=None, line=None):
    print("A warning occurred:")
    print(message)
    print("Do you wish to continue?")

    while True:
        response = input("y/n: ").lower()
        if response not in {"y", "n"}:
            print("Not understood.")
        else:
            break

    if response == "n":
        raise category(message)


# Check if the temporary directory exists, if not, create it
def check_tmp_dir_exists(tmpdir):
    if tmpdir is None:
        tmpdir = tempfile.TemporaryDirectory(dir=os.getcwd())
    else:
        if not os.path.exists(tmpdir):
            _error(f"Temporary directory {tmpdir} does not exist")
            exit(1)
        tmpdir = tempfile.TemporaryDirectory(dir=os.path.abspath(tmpdir))
        # Check if tmpdir has more than 107 characters
        if len(tmpdir.name) > 107:
            _error(f"Temporary directory {tmpdir.name} has more than 107 characters")
            exit(1)
    return tmpdir


def is_debug():
    return bf_logging.should_log(bf_logging.LogLevel.DEBUG)


def _resolve_verbosity_level(
    verbose_count: int, verbose_level: Optional[str]
) -> bf_logging.LogLevel:
    if verbose_level:
        return _VERBOSITY_NAME_TO_LEVEL.get(verbose_level, bf_logging.LogLevel.SUMMARY)

    if verbose_count in _VERBOSITY_COUNT_TO_LEVEL:
        return _VERBOSITY_COUNT_TO_LEVEL[verbose_count]

    if verbose_count <= 0:
        return bf_logging.LogLevel.SUMMARY
    return bf_logging.LogLevel.TRACE


# def refine_chunks(chunks, input_dict, target_weight):
#     # Calculate the current weight of each chunk
#     chunk_weights = [sum(input_dict[key] for key in chunk) for chunk in chunks]

#     while max(chunk_weights) - min(chunk_weights) > target_weight:
#         # Find the chunk with the maximum weight
#         src_chunk = chunk_weights.index(max(chunk_weights))

#         # Find the chunk with the minimum weight
#         dest_chunk = chunk_weights.index(min(chunk_weights))

#         # Find the key in the source chunk with the maximum weight
#         max_weight_key = max(chunks[src_chunk], key=lambda key: input_dict[key])

#         # Move the key with maximum weight from the source to the destination chunk
#         chunks[src_chunk].remove(max_weight_key)
#         chunks[dest_chunk].append(max_weight_key)

#         # Recalculate chunk weights
#         chunk_weights = [sum(input_dict[key] for key in chunk) for chunk in chunks]

#     # Remove empty chunks
#     chunks = [chunk for chunk in chunks if chunk]

#     return chunks


# def sort_keys_by_approx_weight(
#     input_dict, scale=1, num_cores=1, refinement_steps=10, verbose=False
# ):
#     if scale == 0:
#         raise ValueError("Scale cannot be zero.")

#     # Calculate the target weight for each chunk
#     target_weight = scale * int(max(input_dict.values()))

#     # Determine the initial number of chunks based on the number of cores
#     # num_chunks = num_cores * scale
#     num_chunks = (((sum(input_dict.values()) // target_weight)) // num_cores) + 1
#     if num_chunks < num_cores:
#         num_chunks = num_cores
#     # Sort keys by their weights in descending order
#     sorted_keys = sorted(input_dict, key=lambda k: input_dict[k], reverse=True)

#     # Initialize chunks
#     chunks = [[] for _ in range(num_chunks)]
#     total_weight = [0] * num_chunks

#     # Create a progress bar
#     progress_bar = tqdm.tqdm(
#         total=len(sorted_keys),
#         desc="Distributing keys",
#         unit="k",
#         unit_scale=True,
#         unit_divisor=1000,
#         disable=False,
#         leave=False,
#         ncols=80,
#     )

#     # Distribute keys into chunks with weights close to target_weight
#     for key in sorted_keys:
#         # Find the chunk with the least total weight
#         min_chunk_index = min(range(num_chunks), key=lambda i: total_weight[i])

#         # If adding the key doesn't exceed the target weight, add it to the chunk
#         if total_weight[min_chunk_index] + input_dict[key] <= target_weight:
#             chunks[min_chunk_index].append(key)
#             total_weight[min_chunk_index] += input_dict[key]
#         else:
#             # Find the chunk with the weight closest to the target_weight
#             closest_chunk_index = min(
#                 range(num_chunks),
#                 key=lambda i: abs(total_weight[i] + input_dict[key] - target_weight),
#             )
#             chunks[closest_chunk_index].append(key)
#             total_weight[closest_chunk_index] += input_dict[key]

#     # Close the progress bar
#     progress_bar.close()

#     # Initial balance
#     initial_balance = max(len(chunk) for chunk in chunks) - min(
#         len(chunk) for chunk in chunks
#     )

#     # Refinement step
#     for _ in range(refinement_steps):
#         chunks = refine_chunks(chunks, input_dict, target_weight)

#         # Check for improvement in balance
#         current_balance = max(len(chunk) for chunk in chunks) - min(
#             len(chunk) for chunk in chunks
#         )
#         if current_balance >= initial_balance:
#             break  # No improvement, exit the loop

#         # Update initial balance for the next iteration
#         initial_balance = current_balance

#     # Print the min, max, and average weight of each chunk
#     if verbose:
#         for i, chunk in enumerate(chunks, 1):
#             chunk_weights = [input_dict[key] for key in chunk]
#             min_weight = min(chunk_weights)
#             max_weight = max(chunk_weights)
#             avg_weight = sum(chunk_weights) / len(chunk_weights)
#             print(
#                 f"Chunk {i}: Total = {sum(chunk_weights)}, Min Weight = {min_weight}, Max Weight = {max_weight}, Average Weight = {avg_weight}"
#             )

#     return chunks


def refine_chunks(chunks, input_dict, chunk_weights, target_weight):
    while max(chunk_weights) - min(chunk_weights) > target_weight:
        src_idx = chunk_weights.index(max(chunk_weights))
        dest_idx = chunk_weights.index(min(chunk_weights))
        max_weight_key = max(chunks[src_idx], key=lambda key: input_dict[key])

        # Update the weights directly instead of recalculating
        max_key_weight = input_dict[max_weight_key]
        chunk_weights[src_idx] -= max_key_weight
        chunk_weights[dest_idx] += max_key_weight

        # Move the key
        chunks[src_idx].remove(max_weight_key)
        chunks[dest_idx].append(max_weight_key)

    # Remove empty chunks (if any)
    return [chunk for chunk in chunks if chunk]


def sort_keys_by_approx_weight(
    input_dict,
    ref_positions=None,  # Make positions optional
    scale=1,
    num_cores=1,
    refinement_steps=10,
    max_entries_per_chunk=None,
    mode="weight",
    num_entries=1000,
    verbose=False,
):
    if scale == 0:
        raise ValueError("Scale cannot be zero.")

    num_cores = int(scale * num_cores)
    total_weight = sum(input_dict.values())
    if (
        max_entries_per_chunk is not None
        and max_entries_per_chunk > 0
        and max_entries_per_chunk < total_weight
    ):
        target_weight = max_entries_per_chunk
    else:
        target_weight = scale * max(input_dict.values())

    if mode == "weight":
        num_chunks = max(num_cores, int((total_weight // target_weight) + 1))
    else:  # mode == "entries"
        total_entries = len(input_dict)
        num_chunks = max(num_cores, int((total_entries // num_entries) + 1))
        num_chunks = min(num_chunks, 100)  # Cap at 100 chunks

    # Sort items by position and weight if available, else by weight
    if ref_positions is not None:
        sorted_items = sorted(
            ((k, v, ref_positions.get(k, float("inf"))) for k, v in input_dict.items()),
            key=lambda x: (x[2], -x[1]),
        )
        sorted_items = [(item[0], item[1]) for item in sorted_items]
    else:
        sorted_items = sorted(input_dict.items(), key=lambda x: x[1], reverse=True)

    # Greedy bin-packing: assign each key to the chunk with the lowest current weight
    chunks = [[] for _ in range(num_chunks)]
    chunk_weights = [0] * num_chunks
    for key, weight in sorted_items:
        idx = chunk_weights.index(min(chunk_weights))
        chunks[idx].append(key)
        chunk_weights[idx] += weight

    # Refinement: move heaviest item from heaviest to lightest chunk if it improves balance
    for _ in range(refinement_steps):
        max_idx = chunk_weights.index(max(chunk_weights))
        min_idx = chunk_weights.index(min(chunk_weights))
        if max_idx == min_idx:
            break
        # Find the heaviest item in the heaviest chunk
        if not chunks[max_idx]:
            break
        heaviest_key = max(chunks[max_idx], key=lambda k: input_dict[k])
        heaviest_weight = input_dict[heaviest_key]
        # Only move if it improves balance
        if (
            chunk_weights[max_idx] - heaviest_weight
            < chunk_weights[min_idx] + heaviest_weight
        ):
            break
        chunks[max_idx].remove(heaviest_key)
        chunks[min_idx].append(heaviest_key)
        chunk_weights[max_idx] -= heaviest_weight
        chunk_weights[min_idx] += heaviest_weight

    if verbose:
        print_chunk_stats(chunks, input_dict)
        if ref_positions is not None:
            print("\nBAM position stats:")
            for i, chunk in enumerate(chunks):
                if chunk:
                    positions = [ref_positions.get(ref, 0) for ref in chunk]
                    positions.sort()
                    span = positions[-1] - positions[0] if len(positions) > 1 else 0
                    jumps = sum(
                        abs(positions[j] - positions[j - 1])
                        for j in range(1, len(positions))
                    )
                    avg_jump = jumps / (len(positions) - 1) if len(positions) > 1 else 0
                    print(
                        f"Chunk {i+1}: Position span = {span}, Avg ref-to-ref jump = {avg_jump:.2f}, "
                        f"Refs: {chunk[0]}..{chunk[-1]}"
                    )

    return chunks


def print_chunk_stats(chunks, input_dict):
    for i, chunk in enumerate(chunks, 1):
        weights = [input_dict[key] for key in chunk]
        if weights:
            print(
                f"Chunk {i}: Total = {sum(weights)}, "
                f"Min = {min(weights)}, "
                f"Max = {max(weights)}, "
                f"Avg = {sum(weights) / len(chunk)}, "
                f"Size = {len(chunk)}"
            )


def create_empty_output_files(out_files):
    for key, value in out_files.items():
        if value is not None:
            # check if file exists, if not create it
            if os.path.exists(value):
                continue
            if (
                key == "bam_filtered"
                or key == "bam_reassigned"
                or key == "bam_reassigned"
            ):
                create_empty_bam(value)
            elif (
                key == "bam_filtered_tmp"
                or key == "bam_tmp"
                or key == "bam_tmp_sorted"
                or key == "bam_reassigned_tmp"
                or key == "bam_reassigned_sorted"
            ):
                continue
            else:
                Path(value).touch()


# function that creates an empty bam file
def create_empty_bam(output):
    """
    Create an empty bam file
    """
    header = {"HD": {"VN": "1.0", "SO": "unsorted"}}
    # Create an empty BAM file with the specified header
    with pysam.AlignmentFile(output, "wb", header=header) as outfile:
        pass


def check_values(val, minval, maxval, parser, var):
    try:
        value = float(val)
    except ValueError:
        parser.error(
            f"argument {var}: Invalid value {val}. Value has to be a 'float' between {minval} and {maxval}!"
        )
    value = float(val)
    if value < minval or value > maxval:
        parser.error(
            "argument %s: Invalid value %s. Range has to be between %s and %s!"
            % (
                var,
                value,
                minval,
                maxval,
            )
        )
    return float(value)


def check_values_auto(val, minval, maxval, parser, var):
    if val == "auto" or val is None or val == "None":
        return val
    else:
        # check if float
        try:
            val = float(val)
        except ValueError:
            parser.error(
                f"argument {var}: Invalid value {val}. Value has to be 'auto' or a 'float' between {minval} and {maxval}!"
            )
        return check_values(val, minval, maxval, parser, var)


def is_integer(value):
    try:
        int(value)
        return True
    except ValueError:
        return False


def check_suffix(val, parser, var):
    # Define valid units for the argument
    if var == "--scale":
        units = ["K", "M"]
        multiplier = 1000  # Base-1000 scaling
    else:
        units = ["K", "M", "G"]
        multiplier = 1024  # Base-1024 scaling

    # Extract the value and unit
    if val[-1].upper() in units:
        unit = val[-1].upper()
        value = val[:-1]
    else:
        unit = None
        value = val

    # Validate and convert
    if is_integer(value) and int(value) > 0:
        value = int(value)
        if unit == "K":
            value *= multiplier
        elif unit == "M":
            value *= multiplier**2
        elif unit == "G":
            value *= multiplier**3
        return value
    else:
        parser.error(
            f"argument {var}: Invalid value {val}. "
            f"Must be a positive integer optionally followed by {'/'.join(units)}."
        )


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {
        "gz": (b"\x1f", b"\x8b", b"\x08"),
        "bz2": (b"\x42", b"\x5a", b"\x68"),
        "zip": (b"\x50", b"\x4b", b"\x03", b"\x04"),
    }
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, "rb")
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = "plain"
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == "bz2":
        sys.exit("Error: cannot use bzip2 format - use gzip instead")
        sys.exit("Error: cannot use zip format - use gzip instead")
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == "gz":
        return gzip.open
    else:  # plain text
        return open


# From: https://stackoverflow.com/a/11541450
def is_valid_file(parser, arg, var):
    if not os.path.exists(arg):
        parser.error("argument %s: The file %s does not exist!" % (var, arg))
    else:
        return arg


# From https://stackoverflow.com/a/59617044/15704171
def convert_list_to_str(lst):
    n = len(lst)
    if not n:
        return ""
    if n == 1:
        return lst[0]
    return ", ".join(lst[:-1]) + f" or {lst[-1]}"


lca_ranks = [
    "superkingdom",
    "domain",
    "lineage",
    "kingdom",
    "subkingdom",
    "superphylum",
    "phylum",
    "subphylum",
    "superclass",
    "class",
    "subclass",
    "infraclass",
    "clade",
    "cohort",
    "subcohort",
    "superorder",
    "order",
    "suborder",
    "infraorder",
    "parvorder",
    "superfamily",
    "family",
    "subfamily",
    "tribe",
    "subtribe",
    "infratribe",
    "genus",
    "subgenus",
    "section",
    "series",
    "subseries",
    "subsection",
    "species",
    "species group",
    "species subgroup",
    "subspecies",
    "varietas",
    "morph",
    "subvariety",
    "forma",
    "forma specialis",
    "biotype",
    "genotype",
    "isolate",
    "pathogroup",
    "serogroup",
    "serotype",
    "strain",
]


def check_lca_ranks(val, parser, var):
    value = str(val)
    if value in lca_ranks:
        return value
    else:
        parser.error(
            f"argument {var}: Invalid value {value}. Filter has to be one of {convert_list_to_str(lca_ranks)}"
        )


defaults = {
    "min_read_length": 30,
    "max_read_length": 10000,  # ✓ SYNCED: concrete value instead of 0x7fffffff
    "min_read_count": 1,  # ✓ SYNCED: changed from 3 to 1
    "min_expected_breadth_ratio": 0,
    "min_norm_entropy": 0,
    "min_norm_gini": 1.0,
    "min_avg_read_ani": 90.0,
    "min_read_ani": 0.0,
    "min_breadth": 0,
    "min_coverage_evenness": 0,
    "min_coeff_var": float("inf"),
    "min_coverage_mean": 0,
    "prefix": None,
    "sort_memory": "1G",
    "reference_lengths": None,
    "scale": 1_000_000,
    "chunk_size": None,
    "coverage_plots": None,
    "stats": None,
    "stats_filtered": None,
    "bam_filtered": None,
    "bam_reassigned": None,
    "knee_plot": None,
    "read_length_freqs": None,
    "read_hits_count": None,
    "tmp_dir": None,
    "max_em_iterations": 50,
    "em_tolerance": 1e-6,
    "min_probability": 1e-6,
    "prob_fraction": 0,
    "prior_weight": 0.01,
    "use_squarem_acceleration": True,
    "enable_globalization": True,
    "squarem_start_iter": 2,
    "backtrack_factor": 0.5,
    "max_backtrack_steps": 5,
    "steplength_scheme": 3,
    "calculate_pmd": True,
    "library_type": "ds",
    "rank_lca": "genus",
    "lca_summary": None,
    "reference_stats_tsv": None,
    "taxonomy_db": None,
    # Graph construction parameters (cluster-aware filtering always enabled)
    "graph_min_edge_weight": 0,  # Minimum edge weight (shared reads) to keep in graph (0=auto, -1=no filtering)
    # Community detection parameters
    "clustering_algorithm": "leiden",  # "leiden" | "union-find" (simple connected components)
    "community_resolution": 1.0,  # Resolution parameter for community detection
    "community_max_iterations": 10,  # Maximum iterations for convergence
    "community_parallel": False,  # Enable parallel move phase
    "outlier_method": "mad",  # "mad" | "iqr" - Statistical outlier detection method (simplified)
}

help_msg = {
    "bam": "BAM file containing aligned reads",
    "threads": "Number of threads to use",
    "prefix": "Prefix used for the output files",
    "min_read_length": "Minimum read length",
    "max_read_length": "Maximum read length",
    "min_read_count": "Minimum read count",
    "trim_ends": "Exclude n bases at the ends of the reference sequences",
    "trim_min": "Remove coverage that are below this percentile. Used for the Truncated Average Depth (TAD) calculation",
    "trim_max": "Remove coverage that are above this percentile. Used for the Truncated Average Depth (TAD) calculation",
    "min_breadth": "Minimum breadth",
    "min_expected_breadth_ratio": "Minimum expected breadth ratio",
    "min_norm_spatial_entropy": "Minimum normalized spatial entropy (positional distribution uniformity)",
    "min_norm_gini": "Minimum normalized Gini coefficient",
    "min_read_ani": "Minimum read ANI to keep a read",
    "min_avg_read_ani": "Minimum average read ANI",
    "min_coverage_evenness": "Minimum coverage evenness",
    "min_coeff_var": "Minimum coverage evenness calculated as SD/MEAN",
    "min_coverage_mean": "Minimum coverage mean",
    "transform_cov_evenness": "Include those references that fulfill all filtering criteria but the coverage evenness is 0",
    "sort_memory": "Set maximum memory per thread for sorting; suffix K/M/G recognized",
    "scale": "Scale taxonomic abundance by this factor; suffix K/M recognized",
    "read_length_freqs": "Save a JSON file with the read length frequencies mapped to each reference",
    "read_hits_count": "Save a TSV file with the read hits frequencies mapped to each reference",
    "stats": "Save a TSV file with the statistics for each reference",
    "stats_filtered": "Save a TSV file with the statistics for each reference after filtering",
    "bam_filtered": "Save a BAM file with the references that passed the filtering criteria",
    "bam_reassigned": "Save a BAM file without multimapping reads",
    "coverage_plots": "Folder where to save genome coverage plots",
    "knee_plot": "Plot knee plot",
    "sort_by_name": "Sort by read names",
    "disable_sort": "Disable sorting of the filtered BAM file",
    "chunk_size": "Chunk size for parallel processing",
    "tmp_dir": "Temporary directory",
    "help": "Help message",
    "verbose": "Increase logging verbosity (-v for info, -vv for debug, -vvv for trace) or use --verbose LEVEL",
    "reference_lengths": "File with references lengths",
    "low_memory": "Activate the low memory mode",
    "reassign": "Run an EM algorithm to reassign reads to references",
    "reassign_method": "Method for the EM algorithm",
    "reassign_iters": "Number of iterations for the EM algorithm",
    "rank_lca": "Calculate LCA for each read and estimate abundances",
    "names": "Names dmp file from taxonomy",
    "nodes": "Nodes dmp file from taxonomy",
    "acc2taxid": "acc2taxid file from taxonomy",
    "taxonomy_db": "Directory with Parquet taxonomy database (created via filterBAM build-taxonomy, must include accession_map.parquet)",
    "rank_lca": "Rank to use for LCA calculation",
    "lca_summary": "Save a TSV file with the LCA summary",
    "lca_missing": "Save a TSV file with references with missing taxids",
    "lca_stats_flag": "Augment the LCA summary with full per-taxid quality metrics (same output file, more columns).",
    "version": "Print program version",
    # ✓ SYNCED: Updated EM algorithm help messages
    "max_em_iterations": "Maximum number of EM iterations",
    "em_tolerance": "EM convergence tolerance (||F(θ)-θ|| <= ε)",
    "min_probability": "Minimum probability threshold for keeping alignments",
    "prob_fraction": "Relative probability threshold for post-EM filtering",
    "prior_weight": "Prior weight (regularization) for EM algorithm",
    "use_squarem_acceleration": "Enable SQUAREM acceleration for EM algorithm",
    # ✓ NEW: PAPER-SPECIFIC SQUAREM help messages
    "enable_globalization": "Enable gSQUAREM with likelihood monotonicity (backtracking)",
    "squarem_start_iter": "EM iteration to start SQUAREM acceleration",
    "backtrack_factor": "Factor for backtracking α toward -1 (0.1-0.9)",
    "max_backtrack_steps": "Maximum number of backtracking iterations",
    "steplength_scheme": "SQUAREM steplength scheme: 1=S1, 2=S2, 3=S3 (recommended)",
    "output_bam": "Output BAM file with reassigned reads",
    "disable_pmd": "Disable PMD (Post-Mortem Damage) calculation (PMD enabled by default)",
    "library_type": "Library type for PMD scoring: 'ds' (double-stranded) or 'ss' (single-stranded).",
    "single_stranded": "Shortcut for --library-type ss (single-stranded libraries).",
    "reference_stats_tsv": "Save per-reference statistics (TSV). Includes clustering/community columns when clustering is enabled.",
    # Graph construction parameters (cluster-aware filtering always enabled)
    "graph_min_edge_weight": "Minimum edge weight (shared reads between references) to keep edges in graph. 0=auto (default), -1=no filtering, >0=use value. Lower values keep more edges (more connected), higher values prune weak connections (fewer components).",
    # Community clustering help messages
    "clustering_algorithm": "Clustering algorithm for community detection: 'leiden' (high-quality modularity optimization, default) or 'union-find' (fast, simple connected components)",
    "community_resolution": "Resolution parameter for community detection (Leiden γ). Lower=larger communities, higher=smaller communities.",
    "community_max_iterations": "Maximum iterations for the community detection refinement loop (default: 10).",
    "community_parallel": "Enable parallel processing in the community move phase for speed (experimental, Leiden-only).",
    "outlier_method": "Outlier test for clustering coefficients: 'mad' (Median Absolute Deviation, default) or 'iqr' (Interquartile Range). These rules flag structural outliers; multi-metric checks live in the tiered filtering pipeline.",
}

from difflib import get_close_matches, SequenceMatcher


class SubcommandHelpFormatter(argparse.ArgumentParser):
    def _match_score(self, a, b):
        # Remove leading dashes for comparison
        a = a.lstrip("-")
        b = b.lstrip("-")
        return SequenceMatcher(None, a, b).ratio()

    def _get_close_matches(self, arg, possibilities, n=3, cutoff=0.65):
        # Remove leading dashes from the argument for comparison
        clean_arg = arg.lstrip("-")

        # Calculate match scores for all possibilities
        matches = []
        for p in possibilities:
            score = self._match_score(arg, p)
            if score > cutoff:
                matches.append((p, score))

        # Sort by match score and take top N
        matches.sort(key=lambda x: x[1], reverse=True)
        best_matches = [m[0] for m in matches[:n]]

        if best_matches:
            if len(best_matches) == 1:
                return f"\nDid you mean: {best_matches[0]}?"
            else:
                return f"\nDid you mean one of these: {', '.join(best_matches)}?"
        return ""

    def _get_available_commands(self):
        commands = []
        for action in self._actions:
            if isinstance(action, argparse._SubParsersAction):
                commands.extend(action.choices.keys())
        return commands

    def _get_available_arguments(self, parser=None):
        if parser is None:
            parser = self

        arguments = []
        for action in parser._actions:
            if action.option_strings:
                arguments.extend(action.option_strings)
        return sorted(arguments)

    def _print_message(self, message, file=None):
        if message:
            if file is None:
                file = sys.stderr if "error:" in message else sys.stdout
            file.write(message)

    def _clean_error_message(self, message):
        if "unrecognized arguments:" in message:
            parts = message.split(":")
            arg_part = parts[1].strip()
            flag = arg_part.split()[0]
            return f"unrecognized arguments: {flag}"
        return message

    def error(self, message):
        commands = self._get_available_commands()

        try:
            clean_message = self._clean_error_message(message)
            error_msg = f"error: {clean_message}"
            suggestion = ""

            if len(sys.argv) > 1:
                if "invalid choice" in message and sys.argv[1] not in commands:
                    suggestion = self._get_close_matches(
                        sys.argv[1], commands, n=3, cutoff=0.65
                    )
                elif sys.argv[1] in commands:
                    for action in self._actions:
                        if isinstance(action, argparse._SubParsersAction):
                            subparser = action.choices[sys.argv[1]]
                            available_args = self._get_available_arguments(subparser)

                            for arg in sys.argv[2:]:
                                if arg.startswith("-"):
                                    flag = arg.split("=")[0]
                                    if flag not in available_args:
                                        suggestion = self._get_close_matches(
                                            flag, available_args, n=3, cutoff=0.65
                                        )
                                        break

            sys.stderr.write(f"{error_msg}{suggestion}\n")
            sys.exit(2)

        except Exception:
            sys.stderr.write(f"error: {message}\n")
            sys.exit(2)


def get_arguments(argv=None):
    # Create the base parent parser for common arguments
    parent_parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)

    # Add verbosity controls to parent parser
    parent_parser.add_argument(
        "-v",
        dest="verbose_count",
        action="count",
        default=0,
        help=help_msg["verbose"],
    )
    parent_parser.add_argument(
        "--verbose",
        dest="verbose_level",
        choices=tuple(_VERBOSITY_NAME_TO_LEVEL.keys()),
        metavar="LEVEL",
        help="Explicit verbosity level (quiet, summary, info, debug, trace)",
    )

    # Create the main parser
    parser = SubcommandHelpFormatter(
        description="A simple tool to calculate metrics from a BAM file and filter with uneven coverage.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=False,
    )

    # Add version to main parser
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + __version__,
        help=help_msg["version"],
    )

    # Create the shared parser for common required arguments
    common_required = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    required = common_required.add_argument_group("required arguments")
    required.add_argument(
        "--bam",
        required=False,  # Made optional for commands like --list-columns
        dest="bam",
        type=lambda x: is_valid_file(parser, x, "bam"),
        help=help_msg["bam"],
    )

    # Create the shared parser for common optional arguments
    common_optional = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    optional = common_optional.add_argument_group("optional arguments")
    optional.add_argument(
        "-p",
        "--prefix",
        type=str,
        default=defaults["prefix"],
        metavar="STR",
        dest="prefix",
        help=help_msg["prefix"],
    )
    optional.add_argument(
        "-r",
        "--reference-lengths",
        type=lambda x: is_valid_file(parser, x, "reference_lengths"),
        metavar="FILE",
        default=defaults["reference_lengths"],
        dest="reference_lengths_tsv",
        help=help_msg["reference_lengths"],
    )
    optional.add_argument(
        "-t",
        "--threads",
        type=lambda x: int(
            check_values(x, minval=1, maxval=1000, parser=parser, var="--threads")
        ),
        dest="threads",
        metavar="INT",
        default=1,
        help=help_msg["threads"],
    )

    # Create subparsers
    sub_parsers = parser.add_subparsers(
        help="positional arguments",
        dest="action",
    )

    # Create the parser for the filter command with all parent parsers
    parser_filter = sub_parsers.add_parser(
        "filter",
        help="Filter references based on coverage and other metrics",
        parents=[parent_parser, common_required, common_optional],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=False,
    )

    # Create the parser for the reassign command with all parent parsers
    parser_reassign = sub_parsers.add_parser(
        "reassign",
        help="Reassign reads to references using an EM algorithm",
        parents=[parent_parser, common_required, common_optional],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=False,
    )

    # Create the parser for the lca command with all parent parsers
    parser_lca = sub_parsers.add_parser(
        "lca",
        help="Calculate LCA for each read and estimate abundances",
        parents=[parent_parser, common_required, common_optional],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=False,
    )

    # Create the parser for the to-parquet command with all parent parsers
    parser_to_parquet = sub_parsers.add_parser(
        "to-parquet",
        help="Convert BAM to partitioned Parquet format for fast DuckDB queries",
        parents=[parent_parser, common_required, common_optional],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=False,
    )

    # Reassign workflow: organize options by stage
    reassign_io_args = parser_reassign.add_argument_group("Input / Output")
    reassign_pmd_args = parser_reassign.add_argument_group("PMD & Library Type")
    reassign_filter_args = parser_reassign.add_argument_group("Read Filtering")
    reassign_em_args = parser_reassign.add_argument_group("EM Algorithm")
    reassign_squarem_args = parser_reassign.add_argument_group("SQUAREM Acceleration")
    reassign_graph_args = parser_reassign.add_argument_group("Graph Construction")
    reassign_clustering_args = parser_reassign.add_argument_group(
        "Clustering & Outlier Detection"
    )
    reassign_taxonomy_args = parser_reassign.add_argument_group("Taxonomy Integration")
    reassign_export_args = parser_reassign.add_argument_group("Reporting & Export")

    filter_required_args = parser_filter.add_argument_group("Filter required arguments")
    # filter_optional_args = parser_filter.add_argument_group("Filter optional arguments")

    # For LCA, add --taxonomy-db to the existing "required arguments" group from parent
    # Find the existing required arguments group
    lca_required_args = None
    for group in parser_lca._action_groups:
        if group.title == "required arguments":
            lca_required_args = group
            break
    if lca_required_args is None:
        lca_required_args = parser_lca.add_argument_group("required arguments")

    lca_optional_args = parser_lca.add_argument_group("LCA optional arguments")

    # To-parquet argument groups
    parquet_required_args = parser_to_parquet.add_argument_group("required arguments")
    parquet_output_args = parser_to_parquet.add_argument_group("output arguments")
    parquet_format_args = parser_to_parquet.add_argument_group("format options")
    parquet_filter_args = parser_to_parquet.add_argument_group("filtering arguments")

    # add subparser for filtering options:
    # reassign_args = parser.add_argument_group("reassign arguments")
    filtering_filt_args = parser_filter.add_argument_group("filtering arguments")
    # lca_args = parser.add_argument_group("lca arguments")
    misc_filter_args = parser_filter.add_argument_group("miscellaneous arguments")
    out_filter_args = parser_filter.add_argument_group("output arguments")

    reassign_em_args.add_argument(
        "-i",
        "--max-em-iterations",  # ✓ SYNCED: was "--iters"
        type=lambda x: int(
            check_values(
                x, minval=1, maxval=1000, parser=parser, var="--max-em-iterations"
            )
        ),
        metavar="INT",
        default=defaults["max_em_iterations"],  # ✓ SYNCED: new parameter name
        dest="max_em_iterations",  # ✓ SYNCED: new dest name
        help=help_msg["max_em_iterations"],
    )

    reassign_em_args.add_argument(
        "--em-tolerance",
        type=lambda x: float(
            check_values(
                x, minval=1e-12, maxval=1.0, parser=parser, var="--em-tolerance"
            )
        ),
        metavar="FLOAT",
        default=defaults["em_tolerance"],
        dest="em_tolerance",
        help=help_msg["em_tolerance"],
    )

    reassign_em_args.add_argument(
        "--min-probability",  # ✓ SYNCED: was "--min-prob"
        type=lambda x: float(
            check_values(
                x, minval=1e-12, maxval=1.0, parser=parser, var="--min-probability"
            )
        ),
        default=defaults["min_probability"],  # ✓ SYNCED: new parameter name
        metavar="FLOAT",
        dest="min_probability",  # ✓ SYNCED: new dest name
        help=help_msg["min_probability"],
    )

    reassign_em_args.add_argument(
        "--prob-fraction",
        type=lambda x: float(
            check_values(x, minval=0, maxval=1, parser=parser, var="--prob-fraction")
        ),
        default=defaults["prob_fraction"],
        metavar="FLOAT",
        dest="prob_fraction",
        help=help_msg["prob_fraction"],
    )

    reassign_em_args.add_argument(
        "--prior-weight",
        type=lambda x: float(
            check_values(
                x, minval=1e-12, maxval=1.0, parser=parser, var="--prior-weight"
            )
        ),
        metavar="FLOAT",
        default=defaults["prior_weight"],
        dest="prior_weight",
        help=help_msg["prior_weight"],
    )

    # ✓ SYNCED: Read Filtering Parameters
    reassign_filter_args.add_argument(
        "-A",
        "--min-read-ani",
        type=lambda x: float(
            check_values(x, minval=0, maxval=100, parser=parser, var="--min-read-ani")
        ),
        metavar="FLOAT",
        default=defaults["min_read_ani"],
        dest="min_read_ani",
        help=help_msg["min_read_ani"],
    )

    reassign_filter_args.add_argument(
        "-l",
        "--min-read-length",
        type=lambda x: int(
            check_values(
                x, minval=1, maxval=100000, parser=parser, var="--min-read-length"
            )
        ),
        default=defaults["min_read_length"],
        metavar="INT",
        dest="min_read_length",
        help=help_msg["min_read_length"],
    )

    reassign_filter_args.add_argument(
        "-L",
        "--max-read-length",
        type=lambda x: int(
            check_values(
                x,
                minval=1,
                maxval=100000,
                parser=parser,
                var="--max-read-length",  # ✓ SYNCED: concrete max
            )
        ),
        default=defaults["max_read_length"],  # ✓ SYNCED: now 10000
        metavar="INT",
        dest="max_read_length",
        help=help_msg["max_read_length"],
    )

    reassign_filter_args.add_argument(
        "-n",
        "--min-read-count",
        type=lambda x: int(
            check_values(
                x, minval=1, maxval=float("inf"), parser=parser, var="--min-read-count"
            )
        ),
        default=defaults["min_read_count"],  # ✓ SYNCED: now 1
        metavar="INT",
        dest="min_read_count",
        help=help_msg["min_read_count"],
    )

    reassign_squarem_args.add_argument(
        "--disable-squarem",
        dest="use_squarem_acceleration",
        action="store_false",
        default=True,  # Enabled by default
        help="Disable SQUAREM acceleration (use standard EM instead)",
    )

    reassign_squarem_args.add_argument(
        "--disable-globalization",
        dest="enable_globalization",
        action="store_false",
        default=True,  # Enabled by default
        help="Disable gSQUAREM globalization (use non-monotone SQUAREM)",
    )

    reassign_squarem_args.add_argument(
        "--squarem-start-iter",  # ✓ NEW: when to start SQUAREM
        type=lambda x: int(
            check_values(
                x, minval=1, maxval=50, parser=parser, var="--squarem-start-iter"
            )
        ),
        default=defaults["squarem_start_iter"],
        metavar="INT",
        dest="squarem_start_iter",
        help=help_msg["squarem_start_iter"],
    )

    reassign_squarem_args.add_argument(
        "--backtrack-factor",  # ✓ NEW: backtracking control
        type=lambda x: float(
            check_values(
                x, minval=0.1, maxval=0.9, parser=parser, var="--backtrack-factor"
            )
        ),
        default=defaults["backtrack_factor"],
        metavar="FLOAT",
        dest="backtrack_factor",
        help=help_msg["backtrack_factor"],
    )

    reassign_squarem_args.add_argument(
        "--max-backtrack-steps",  # ✓ NEW: backtracking limit
        type=lambda x: int(
            check_values(
                x, minval=1, maxval=20, parser=parser, var="--max-backtrack-steps"
            )
        ),
        default=defaults["max_backtrack_steps"],
        metavar="INT",
        dest="max_backtrack_steps",
        help=help_msg["max_backtrack_steps"],
    )

    reassign_squarem_args.add_argument(
        "--steplength-scheme",  # ✓ NEW: S1/S2/S3 selection
        type=lambda x: int(
            check_values(
                x, minval=1, maxval=3, parser=parser, var="--steplength-scheme"
            )
        ),
        default=defaults["steplength_scheme"],
        metavar="INT",
        dest="steplength_scheme",
        help=help_msg["steplength_scheme"],
    )

    # ✓ SYNCED: Output Parameters
    reassign_io_args.add_argument(
        "-o",
        "--output-bam",  # ✓ SYNCED: was "--out-bam"
        dest="output_bam",  # ✓ SYNCED: new dest name
        default=defaults["bam_reassigned"],
        metavar="FILE",
        type=str,
        nargs="?",
        const="",
        help=help_msg["output_bam"],
    )
    reassign_pmd_args.add_argument(
        "--library-type",
        dest="library_type",
        choices=("ds", "ss"),
        default=defaults["library_type"],
        metavar="{ds,ss}",
        help=help_msg["library_type"],
    )
    reassign_pmd_args.add_argument(
        "--single-stranded",
        dest="library_type",
        action="store_const",
        const="ss",
        default=argparse.SUPPRESS,
        help=help_msg["single_stranded"],
    )
    reassign_pmd_args.add_argument(
        "--double-stranded",
        dest="library_type",
        action="store_const",
        const="ds",
        default=argparse.SUPPRESS,
        help="Shortcut for --library-type ds (double-stranded libraries).",
    )
    reassign_pmd_args.add_argument(
        "--disable-pmd",
        dest="disable_pmd",
        action="store_true",
        default=False,  # PMD enabled by default
        help=help_msg["disable_pmd"],
    )
    reassign_export_args.add_argument(
        "-S",
        "--reference-stats",
        dest="reference_stats_tsv",
        default=defaults["reference_stats_tsv"],
        metavar="FILE",
        type=str,
        nargs="?",
        const="",
        help=help_msg["reference_stats_tsv"],
    )
    # Graph construction arguments (cluster-aware filtering always enabled)
    reassign_graph_args.add_argument(
        "--graph-min-edge-weight",
        dest="graph_min_edge_weight",
        type=lambda x: (
            int(x)
            if int(x) >= -1
            else parser.error(
                f"argument --graph-min-edge-weight: Invalid value {x}. Must be -1 (no filtering), 0 (auto), or positive integer."
            )
        ),
        default=defaults["graph_min_edge_weight"],
        metavar="INT",
        help=help_msg["graph_min_edge_weight"],
    )

    # Clustering flags
    reassign_clustering_args.add_argument(
        "--clustering",
        dest="clustering",
        action="store_true",
        default=False,
        help="Enable clustering / community detection (requires --reference-stats)",
    )
    # Keep community detection tuning parameters (resolution / max iterations) but do not expose
    # clustering algorithm selection or parallel flag via CLI to simplify the interface.
    reassign_clustering_args.add_argument(
        "--community-resolution",
        dest="community_resolution",
        type=float,
        default=defaults["community_resolution"],
        metavar="FLOAT",
        help=help_msg["community_resolution"],
    )
    reassign_clustering_args.add_argument(
        "--community-max-iterations",
        dest="community_max_iterations",
        type=int,
        default=defaults["community_max_iterations"],
        metavar="INT",
        help=help_msg["community_max_iterations"],
    )
    reassign_clustering_args.add_argument(
        "--outlier-method",
        dest="outlier_method",
        type=str,
        choices=["mad", "iqr"],
        default=defaults["outlier_method"],
        metavar="METHOD",
        help=help_msg["outlier_method"],
    )

    reassign_export_args.add_argument(
        "--graph-export",
        dest="graph_export",
        type=str,
        default=None,
        metavar="PATH",
        help="Export graph to GraphML format for visualization (Cytoscape, igraph, etc.). "
        "Includes node attributes (reference name, community, CC, read counts) and edge weights.",
    )

    # Taxonomy-aware graph analysis
    reassign_taxonomy_args.add_argument(
        "--taxonomy-db",
        dest="taxonomy_db",
        type=str,
        default=None,
        metavar="DIR",
        help="Path to taxonomy database directory (Parquet format). "
        "Should contain nodes.parquet, metadata.parquet, and optionally accession_map.parquet. "
        "If accession_map.parquet exists, taxonomy-aware graph analysis will be enabled.",
    )

    reassign_taxonomy_args.add_argument(
        "--taxonomy-min-rank",
        dest="taxonomy_min_rank",
        type=int,
        default=6,
        metavar="RANK_ID",
        help="Minimum taxonomic rank ID for mismatch detection. "
        "Rank IDs: 2=species, 6=genus, 8=family, 13=order, 17=class, 20=phylum, 24=superkingdom. "
        "Default: 6 (genus level).",
    )

    reassign_taxonomy_args.add_argument(
        "--taxonomy-cross-domain-threshold",
        dest="taxonomy_cross_domain_threshold",
        type=float,
        default=0.10,
        metavar="FRAC",
        help="Fraction of cross-domain neighbors to flag reference as contamination. "
        "Default: 0.10 (flag if >10%% of neighbors cross domains).",
    )

    reassign_taxonomy_args.add_argument(
        "--taxonomy-kingdom-threshold",
        dest="taxonomy_kingdom_threshold",
        type=float,
        default=0.25,
        metavar="FRAC",
        help="Fraction of kingdom-mismatch neighbors to flag reference. "
        "Default: 0.25 (flag if >25%% of neighbors cross kingdoms).",
    )

    reassign_taxonomy_args.add_argument(
        "--taxonomy-genus-threshold",
        dest="taxonomy_genus_threshold",
        type=float,
        default=0.50,
        metavar="FRAC",
        help="Fraction of genus-level mismatch neighbors to flag reference. "
        "Default: 0.50 (flag if >50%% of neighbors differ at genus level).",
    )

    # Taxonomy-informed filtering options (combine graph + taxonomy for automated filtering)
    reassign_taxonomy_args.add_argument(
        "--taxonomy-filter",
        dest="taxonomy_filter_enabled",
        action="store_true",
        help="Enable taxonomy-informed filtering (combines graph topology + taxonomy for better filtering decisions). "
        "Requires --taxonomy-db and --clustering.",
    )

    reassign_taxonomy_args.add_argument(
        "--taxonomy-strict-filter",
        dest="taxonomy_strict_filter",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Enable strict taxonomy filtering for --remove-cross-domain-references mode. "
        "Automatically removes references with cross-domain flags AND high connectivity. "
        "Default: enabled when --taxonomy-filter is set.",
    )

    reassign_taxonomy_args.add_argument(
        "--taxonomy-min-connections",
        dest="taxonomy_strict_min_connections",
        type=int,
        default=5,
        metavar="N",
        help="Minimum graph connections required for reference removal in --remove-cross-domain-references mode. "
        "Default: 5 (only remove cross-domain refs with >= 5 neighbors).",
    )

    reassign_taxonomy_args.add_argument(
        "--taxonomy-weighted-outlier",
        dest="taxonomy_weighted_outlier",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Weight anomaly scores by taxonomy flags (increases sensitivity to taxonomic incongruence). "
        "Default: enabled when --taxonomy-filter is set.",
    )

    reassign_taxonomy_args.add_argument(
        "--taxonomy-anomaly-weight",
        dest="taxonomy_anomaly_weight",
        type=float,
        default=2.0,
        metavar="WEIGHT",
        help="Weight multiplier for taxonomy-flagged references in outlier detection. "
        "Default: 2.0 (2x more suspicious if taxonomy flag is set).",
    )

    reassign_taxonomy_args.add_argument(
        "--taxonomy-second-chance",
        dest="taxonomy_second_chance",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Give graph-flagged references with normal taxonomy a second chance (reduces false positives). "
        "Default: enabled when --taxonomy-filter is set.",
    )

    reassign_taxonomy_args.add_argument(
        "--taxonomy-second-chance-cc",
        dest="taxonomy_second_chance_cc",
        type=float,
        default=0.3,
        metavar="CC",
        help="Minimum clustering coefficient for second-chance validation. "
        "Default: 0.3 (restore removed refs with CC >= 0.3 and normal taxonomy).",
    )

    reassign_taxonomy_args.add_argument(
        "--no-cross-domain-removal",
        dest="no_cross_domain_removal",
        action="store_true",
        help="Disable automatic cross-domain removal when using --taxonomy-filter. "
        "By default, combined cross-domain removal (references + alignments) is enabled automatically "
        "when --taxonomy-filter and --taxonomy-db are used. Use this flag to disable that behavior.",
    )

    reassign_taxonomy_args.add_argument(
        "--remove-cross-domain-alignments",
        dest="remove_cross_domain_alignments",
        action="store_true",
        default=None,
        help="Remove alignments that connect references from different domains (Bacteria/Archaea/Eukaryota). "
        "This preserves references but removes the spurious cross-domain read mappings. "
        "NOTE: This is now enabled by default with --taxonomy-filter. Use --no-cross-domain-removal to disable.",
    )

    reassign_taxonomy_args.add_argument(
        "--remove-cross-domain-references",
        dest="remove_cross_domain_references",
        action="store_true",
        default=None,
        help="Remove entire references that are flagged as cross-domain contamination. "
        "References with high connectivity to other domains are removed along with all their alignments. "
        "NOTE: This is now enabled by default with --taxonomy-filter. Use --no-cross-domain-removal to disable.",
    )

    reassign_taxonomy_args.add_argument(
        "--remove-cross-domain-all",
        dest="remove_cross_domain_all",
        action="store_true",
        help="Apply both alignment-level AND reference-level cross-domain removal (most aggressive). "
        "First removes cross-domain references, then removes remaining cross-domain alignments. "
        "NOTE: This is now the default when using --taxonomy-filter with --taxonomy-db.",
    )

    reassign_taxonomy_args.add_argument(
        "--detect-misannotations",
        dest="detect_misannotations",
        action="store_true",
        default=None,
        help="Detect and flag potential database misannotations by analyzing cross-domain patterns. "
        "References that lose ALL edges after cross-domain removal are flagged with high confidence. "
        "Results exported to TSV with misannotation_flag and confidence_score columns for database curation. "
        "NOTE: This is now enabled by default with --taxonomy-filter. Use --no-cross-domain-removal to disable.",
    )

    misc_filter_args.add_argument(
        "--reference-trim-length",
        type=lambda x: int(
            check_values(
                x, minval=0, maxval=10000, parser=parser, var="---reference-trim-length"
            )
        ),
        dest="trim_ends",
        metavar="INT",
        default=0,
        help=help_msg["trim_ends"],
    )
    misc_filter_args.add_argument(
        "--trim-min",
        type=lambda x: int(
            check_values(x, minval=0, maxval=100, parser=parser, var="--trim-min")
        ),
        dest="trim_min",
        metavar="INT",
        default=10,
        help=help_msg["trim_min"],
    )
    misc_filter_args.add_argument(
        "--trim-max",
        type=lambda x: int(
            check_values(x, minval=0, maxval=100, parser=parser, var="--trim-max")
        ),
        dest="trim_max",
        metavar="INT",
        default=90,
        help=help_msg["trim_max"],
    )
    # Generic column-based filtering - replaces all individual filter arguments
    filtering_filt_args.add_argument(
        "-f",
        "--filter",
        type=str,
        metavar="SPEC",
        dest="filter_spec",
        default=None,
        help=(
            "Generic column-based filtering. Format: 'column:min:max' where "
            "min/max can be empty for unbounded. Multiple filters separated by commas. "
            "Examples: "
            "'read_ani_mean:90:' (>= 90), "
            "'breadth::0.95' (<= 0.95), "
            "'coverage_mean:5:1000' (between 5 and 1000). "
            "Use --list-columns to see all filterable columns."
        ),
    )
    filtering_filt_args.add_argument(
        "--list-columns",
        action="store_true",
        dest="list_columns",
        help="List all available filterable columns with their data types and exit",
    )

    # Keep read filtering arguments (these are read-level filters, not reference-level)
    filtering_filt_args.add_argument(
        "-A",
        "--min-read-ani",
        type=lambda x: float(
            check_values(x, minval=0, maxval=100, parser=parser, var="--min-read-ani")
        ),
        metavar="FLOAT",
        default=defaults["min_read_ani"],
        dest="min_read_ani",
        help="Minimum read ANI to include in statistics (read-level filter, not reference-level)",
    )
    filtering_filt_args.add_argument(
        "-l",
        "--min-read-length",
        type=lambda x: int(
            check_values(
                x, minval=1, maxval=100000, parser=parser, var="--min-read-length"
            )
        ),
        default=defaults["min_read_length"],
        metavar="INT",
        dest="min_read_length",
        help="Minimum read length to include in statistics (read-level filter, not reference-level)",
    )
    filtering_filt_args.add_argument(
        "-L",
        "--max-read-length",
        type=lambda x: int(
            check_values(
                x, minval=1, maxval=float("inf"), parser=parser, var="--max-read-length"
            )
        ),
        default=defaults["max_read_length"],
        metavar="INT",
        dest="max_read_length",
        help="Maximum read length to include in statistics (read-level filter, not reference-level)",
    )

    # Reference-level filter for minimum read count (commonly used, kept as shorthand)
    filtering_filt_args.add_argument(
        "-n",
        "--min-read-count",
        type=lambda x: int(
            check_values(
                x, minval=1, maxval=float("inf"), parser=parser, var="--min-read-count"
            )
        ),
        default=defaults["min_read_count"],
        metavar="INT",
        dest="min_read_count",
        help="Minimum number of reads per reference to pass filter (reference-level filter). "
             "Equivalent to --filter 'read_count:N:' but provided for convenience.",
    )
    misc_filter_args.add_argument(
        "-m",
        "--sort-memory",
        type=lambda x: check_suffix(x, parser=parser, var="--sort-memory"),
        default=defaults["sort_memory"],
        metavar="STR",
        dest="sort_memory",
        help=help_msg["sort_memory"],
    )
    misc_filter_args.add_argument(
        "-N",
        "--sort-by-name",
        dest="sort_by_name",
        action="store_true",
        help=help_msg["sort_by_name"],
    )
    misc_filter_args.add_argument(
        "--disable-sort",
        dest="disable_sort",
        action="store_true",
        help=help_msg["disable_sort"],
    )
    misc_filter_args.add_argument(
        "--scale",
        type=lambda x: check_suffix(x, parser=parser, var="--scale"),
        default=defaults["scale"],
        dest="scale",
        metavar="STR",
        help=help_msg["scale"],
    )
    filter_required_args.add_argument(
        "--stats",
        dest="output",
        default=defaults["stats"],
        type=str,
        metavar="FILE",
        nargs="?",
        const="",
        required=False,  # Made optional for commands like --list-columns
        help=help_msg["stats"],
    )
    out_filter_args.add_argument(
        "--stats-filtered",
        dest="filtered_output",
        default=defaults["stats_filtered"],
        type=str,
        metavar="FILE",
        nargs="?",
        const="",
        help=help_msg["stats_filtered"],
    )
    out_filter_args.add_argument(
        "--bam-filtered",
        dest="filtered_bam",
        default=defaults["bam_filtered"],
        metavar="FILE",
        type=str,
        nargs="?",
        const="",
        help=help_msg["bam_filtered"],
    )
    # parser.add_argument(
    #     "--chunk-size",
    #     type=lambda x: int(
    #         check_values(x, minval=1, maxval=100000, parser=parser, var="--chunk-size")
    #     ),
    #     default=defaults["chunk_size"],
    #     metavar="INT",
    #     dest="chunk_size",
    #     help=help_msg["chunk_size"],
    # )

    lca_required_args.add_argument(
        "--taxonomy-db",
        metavar="DIR",
        type=str,
        dest="taxonomy_db",
        required=True,
        help=help_msg["taxonomy_db"],
    )
    lca_optional_args.add_argument(
        "--lca-rank",
        metavar="STR",
        type=lambda x: str(check_lca_ranks(x, parser=parser, var="--lca-rank")),
        default=defaults["rank_lca"],
        dest="rank_lca",
        help=help_msg["rank_lca"],
    )
    lca_required_args.add_argument(
        "--lca-summary",
        dest="lca_summary",
        metavar="FILE",
        type=str,
        help=help_msg["lca_summary"],
    )
    lca_optional_args.add_argument(
        "--lca-per-read",
        dest="lca_per_read",
        metavar="FILE",
        default=None,
        type=str,
        help="Write per-read LCA assignments (TSV). Use .gz to compress (e.g., per-read.tsv.gz).",
    )
    lca_optional_args.add_argument(
        "--stats",
        dest="lca_stats",
        action="store_true",
        help=help_msg["lca_stats_flag"],
    )

    # Create the parser for the build-taxonomy command
    parser_build_taxonomy = sub_parsers.add_parser(
        "build-taxonomy",
        help="Build Parquet taxonomy database from NCBI dump files",
        parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=False,
    )

    build_taxonomy_required_args = parser_build_taxonomy.add_argument_group(
        "Build taxonomy required arguments"
    )
    build_taxonomy_optional_args = parser_build_taxonomy.add_argument_group(
        "Build taxonomy optional arguments"
    )

    build_taxonomy_required_args.add_argument(
        "--nodes",
        dest="nodes",
        type=str,
        required=True,
        metavar="FILE",
        help="Path to NCBI nodes.dmp file (plain text or .gz compressed)",
    )

    build_taxonomy_required_args.add_argument(
        "--names",
        dest="names",
        type=str,
        required=True,
        metavar="FILE",
        help="Path to NCBI names.dmp file (plain text or .gz compressed)",
    )

    build_taxonomy_required_args.add_argument(
        "--output",
        "-o",
        dest="output",
        type=str,
        required=True,
        metavar="DIR",
        help="Output directory for Parquet taxonomy database files",
    )

    build_taxonomy_required_args.add_argument(
        "--acc2taxid",
        dest="acc2taxid",
        type=str,
        nargs="+",
        required=True,
        metavar="FILE",
        help="Path(s) to NCBI accession2taxid file(s) (plain text or .gz compressed). Can specify multiple files.",
    )

    build_taxonomy_optional_args.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=1,
        metavar="INT",
        help="Number of threads for parallel processing",
    )

    build_taxonomy_optional_args.add_argument(
        "--cache-taxids",
        dest="cache_taxids",
        type=int,
        nargs="+",
        default=None,
        metavar="TAXID",
        help="List of taxids to precompute in LCA cache for O(1) queries. Optional.",
    )

    # To-parquet arguments
    parquet_output_args.add_argument(
        "-o",
        "--output",
        dest="output",
        type=str,
        required=True,
        metavar="DIR",
        help="Output directory for Parquet files (will create partitioned structure)",
    )

    parquet_format_args.add_argument(
        "--num-partitions",
        dest="num_partitions",
        type=int,
        default=-1,
        metavar="INT",
        help="Number of hash-based partitions for alignments (default: auto)",
    )

    parquet_format_args.add_argument(
        "--batch-size",
        dest="batch_size",
        type=int,
        default=-1,
        metavar="INT",
        help="Number of alignments per batch for streaming (default: auto)",
    )

    parquet_format_args.add_argument(
        "--compression",
        dest="compression",
        type=str,
        choices=["snappy", "zstd", "gzip", "none"],
        default="zstd",
        metavar="CODEC",
        help="Parquet compression codec: snappy (fast), zstd (balanced), gzip (max compression) (default: zstd)",
    )

    parquet_format_args.add_argument(
        "--compression-level",
        dest="compression_level",
        type=int,
        default=3,
        metavar="INT",
        help="Compression level for zstd (1-9, default: 3)",
    )

    parquet_format_args.add_argument(
        "--include-read-names",
        dest="include_read_names",
        action="store_true",
        default=True,
        help="Include read names in output (default: enabled; disable with --no-read-names)",
    )

    parquet_format_args.add_argument(
        "--no-read-names",
        dest="include_read_names",
        action="store_false",
        help=argparse.SUPPRESS,
    )

    parquet_format_args.add_argument(
        "--no-sequences",
        dest="include_sequences",
        action="store_false",
        default=True,
        help="Exclude sequences from output (default: sequences included)",
    )

    parquet_format_args.add_argument(
        "--include-sequences",
        dest="include_sequences",
        action="store_true",
        help=argparse.SUPPRESS,
    )

    parquet_format_args.add_argument(
        "--include-sequence-text",
        dest="include_sequence_text",
        action="store_true",
        default=False,
        help="Also store plain-text sequences alongside the packed representation (default: disabled)",
    )

    parquet_format_args.add_argument(
        "--no-sequence-text",
        dest="include_sequence_text",
        action="store_false",
        help=argparse.SUPPRESS,
    )

    parquet_format_args.add_argument(
        "--disable-pmd",
        dest="calculate_pmd",
        action="store_false",
        default=True,
        help="Disable PMD (Post-Mortem Damage) score calculation (enabled by default)",
    )

    parquet_filter_args.add_argument(
        "--min-mapq",
        dest="min_mapq",
        type=int,
        default=0,
        metavar="INT",
        help="Minimum mapping quality (default: 0, no filter)",
    )

    parquet_filter_args.add_argument(
        "--min-read-length",
        dest="min_read_length",
        type=int,
        default=0,
        metavar="INT",
        help="Minimum read length (default: 0)",
    )

    parquet_filter_args.add_argument(
        "--max-read-length",
        dest="max_read_length",
        type=int,
        default=2147483647,
        metavar="INT",
        help="Maximum read length (default: 2,147,483,647; effectively unlimited)",
    )

    parquet_filter_args.add_argument(
        "--min-read-ani",
        dest="min_read_ani",
        type=float,
        default=0.0,
        metavar="FLOAT",
        help="Minimum alignment identity percentage (default: 0.0)",
    )

    if argv is None:
        argv = sys.argv[1:]

    if not argv:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args(argv)

    verbosity_level = _resolve_verbosity_level(
        getattr(args, "verbose_count", 0), getattr(args, "verbose_level", None)
    )
    args.verbosity = verbosity_level
    args.verbose = verbosity_level >= bf_logging.LogLevel.INFO
    args.debug = verbosity_level >= bf_logging.LogLevel.DEBUG
    bf_logging.set_level(verbosity_level)

    # Enforce that clustering requires graph-analysis TSV output
    if getattr(args, "clustering", False) and not getattr(
        args, "reference_stats_tsv", None
    ):
        parser.error("--clustering requires --reference-stats to be set")

    # Print chosen graph edge-weight mode early so users see whether auto/none/value was selected
    try:
        gm = getattr(args, "graph_min_edge_weight", None)
        if gm is not None and args.verbosity >= bf_logging.LogLevel.INFO:
            if gm == -1:
                bf_logging.log("GRAPH", "Mode: none (no filtering). CLI value=%s", gm)
            elif gm == 0:
                bf_logging.log(
                    "GRAPH",
                    "Mode: auto (will choose threshold later). CLI value=%s",
                    gm,
                )
            else:
                bf_logging.log("GRAPH", "Mode: explicit value. CLI value=%s", gm)
    except Exception:
        # Fail silently if args structure unexpected
        pass

    return args


@contextmanager
def suppress_stdout():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(devnull, "w") as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)


def fast_flatten(input_list):
    return list(chain.from_iterable(input_list))


# def concat_df(frames):
#     COLUMN_NAMES = frames[0].columns
#     df_dict = dict.fromkeys(COLUMN_NAMES, [])
#     for col in COLUMN_NAMES:
#         extracted = (frame[col] for frame in frames)
#         # Flatten and save to df_dict
#         df_dict[col] = fast_flatten(extracted)
#     df = pd.DataFrame.from_dict(df_dict)[COLUMN_NAMES]
#     return df


def initializer(init_data):
    global parms
    parms = init_data


def clean_up(keep, temp_dir):
    if keep:
        _info("Cleaning up temporary files")
        shutil.rmtree(temp_dir, ignore_errors=True)


# from https://stackoverflow.com/questions/53751050/python-multiprocessing-understanding-logic-behind-chunksize/54032744#54032744
def calc_chunksize(n_workers, len_iterable, factor=4):
    """Calculate chunksize argument for Pool-methods.

    Resembles source-code within `multiprocessing.pool.Pool._map_async`.
    """
    chunksize, extra = divmod(len_iterable, n_workers * factor)
    if extra:
        chunksize += 1
    return chunksize


# def create_output_files(
#     prefix,
#     bam,
#     stats,
#     stats_filtered,
#     bam_filtered,
#     read_length_freqs,
#     read_hits_count,
#     knee_plot,
#     coverage_plots,
# ):
#     if prefix is None:
#         prefix = bam.replace(".bam", "")

#     out_files = {}
#     if stats is not None:
#         if stats == "":
#             out_files["stats"] = f"{prefix}_stats.tsv.gz"
#         else:
#             out_files["stats"] = stats
#     if stats_filtered is not None:
#         if stats_filtered == "":
#             out_files["stats_filtered"] = f"{prefix}_stats-filtered.tsv.gz"
#         else:
#             out_files["stats_filtered"] = stats_filtered
#     if bam_filtered is not None:
#         if bam_filtered == "":
#             out_files["bam_filtered"] = f"{prefix}.filtered.bam"
#         else:
#             out_files["bam_filtered"] = bam_filtered
#     if read_length_freqs is not None:
#         if read_length_freqs == "":
#             out_files["read_length_freqs"] = f"{prefix}_read-length-freqs.json"
#         else:
#             out_files["read_length_freqs"] = read_length_freqs
#     if read_hits_count is not None:
#         if read_hits_count == "":
#             out_files["read_hits_count"] = f"{prefix}_read-hits-count.tsv.gz"
#         else:
#             out_files["read_hits_count"] = read_hits_count
#     if knee_plot is not None:
#         if knee_plot == "":
#             out_files["knee_plot"] = f"{prefix}_knee-plot.png"
#         else:
#             out_files["knee_plot"] = knee_plot
#     if coverage_plots is not None:
#         if coverage_plots == "":
#             out_files["coverage_plot_dir"] = f"{prefix}_coverage-plots"
#         else:
#             out_files["coverage_plot_dir"] = coverage_plots
#     out_files["bam_filtered_tmp"] = (f"{prefix}.filtered.tmp.bam",)


#     # create output files
#     out_files = {
#         "stats": stats,
#         "stats_filtered": stats_filtered,
#         "bam_filtered_tmp": f"{prefix}.filtered.tmp.bam",
#         "bam_filtered": bam_filtered,
#         "read_length_freqs": read_length_freqs,
#         "read_hits_count": read_hits_count,
#         "knee_plot": knee_plot,
#         "coverage_plot_dir": coverage_plots,
#     }
#     return out_files
def create_output_files(
    bam,
    tmp_dir,
    prefix=None,
    mode=None,
    stats="",
    stats_filtered="",
    bam_reassigned="",
    bam_filtered="",
    read_length_freqs="",
    read_hits_count="",
    knee_plot="",
    coverage_plots="",
    lca_summary="",
):
    if prefix is None:
        prefix = Path(bam).with_suffix("").name

    # Handle temporary directory
    if tmp_dir is None:
        # Create a temporary directory in current working directory
        tmp_dir_path = tempfile.mkdtemp(dir=os.getcwd())
    elif isinstance(tmp_dir, str):
        # User provided a path as string
        if not os.path.exists(tmp_dir):
            _error(f"Temporary directory {tmp_dir} does not exist")
            exit(1)
        tmp_dir_path = os.path.abspath(tmp_dir)
    else:
        # tmp_dir is a TemporaryDirectory object
        tmp_dir_path = tmp_dir.name

    if stats == "" or stats is None:
        stats = f"{prefix}_stats.tsv.gz"
    if stats_filtered == "" or stats_filtered is None:
        stats_filtered = f"{prefix}_stats-filtered.tsv.gz"
    if bam_filtered == "" or bam_filtered is None:
        bam_filtered = f"{prefix}.filtered.bam"
    if bam_reassigned == "" or bam_reassigned is None:
        bam_reassigned = f"{prefix}.reassigned.bam"
    if read_length_freqs == "" or read_length_freqs is None:
        read_length_freqs = f"{prefix}_read-length-freqs.json"
    if read_hits_count == "" or read_hits_count is None:
        read_hits_count = f"{prefix}_read-hits-count.tsv.gz"
    if knee_plot == "" or knee_plot is None:
        knee_plot = f"{prefix}_knee-plot.png"
    if coverage_plots == "" or coverage_plots is None:
        coverage_plots = f"{prefix}_coverage-plots"
    if lca_summary == "" or lca_summary is None:
        lca_summary = f"{prefix}_lca-summary.tsv.gz"

    # create output files
    if mode == "filter":
        out_files = {
            "stats": stats,
            "stats_filtered": stats_filtered,
            "bam_filtered_tmp": f"{tmp_dir_path}/{prefix}.filtered.tmp.bam",
            "bam_filtered": bam_filtered,
            "read_length_freqs": read_length_freqs,
            "read_hits_count": read_hits_count,
            "knee_plot": knee_plot,
            "coverage_plot_dir": coverage_plots,
            "bam_tmp": f"{tmp_dir_path}/{prefix}.tmp.bam",
            "bam_tmp_sorted": f"{tmp_dir_path}/{prefix}.tmp.sorted.bam",
        }
    elif mode == "reassign":
        out_files = {
            "bam_reassigned_tmp": f"{tmp_dir_path}/{prefix}.reassigned.tmp.bam",
            "bam_reassigned_sorted": f"{tmp_dir_path}/{prefix}.reassigned.sorted.bam",
            "bam_reassigned": bam_reassigned,
        }
    elif mode == "lca":
        out_files = {
            "lca_summary": lca_summary,
        }
    else:
        _error("Mode not recognized")
        exit(1)
    out_files["tmp_dir"] = tmp_dir_path
    out_files["sorted_bam"] = f"{tmp_dir_path}/{prefix}.bf-sorted.bam"

    # check that read_length_freqs is a json file
    if read_length_freqs is not None:
        if not read_length_freqs.endswith(".json"):
            _error("--read-length-freqs must be a JSON file")
            exit(1)
    return out_files

    # out_files = {
    #     "stats": stats,
    #     "stats_filtered": stats_filtered,
    #     "bam_filtered_tmp": f"{tmp_dir}/{prefix}.filtered.tmp.bam",
    #     "bam_tmp": f"{tmp_dir}/{prefix}.tmp.bam",
    #     "bam_tmp_sorted": f"{tmp_dir}/{prefix}.tmp.sorted.bam",
    #     "bam_filtered": bam_filtered,
    #     "bam_reassigned_tmp": f"{tmp_dir}/{prefix}.reassigned.tmp.bam",
    #     "bam_reassigned_sorted": f"{tmp_dir}/{prefix}.reassigned.sorted.bam",
    #     "bam_reassigned": bam_reassigned,
    #     "read_length_freqs": read_length_freqs,
    #     "read_hits_count": read_hits_count,
    #     "knee_plot": knee_plot,
    #     "coverage_plot_dir": coverage_plots,
    #     "lca_summary": lca_summary,
    # }


def allocate_threads(total_threads, min_io_processes, max_io_processes):
    """
    Allocates threads between CPU-bound workers and I/O-bound processes based on total threads and
    desired range of I/O processes. Selects the best compromise to maximize CPU-bound workers.

    Parameters:
        total_threads (int): Total number of available threads.
        min_io_processes (int): Minimum desired I/O-bound processes.
        max_io_processes (int): Maximum desired I/O-bound processes.

    Returns:
        tuple: (Number of workers, Number of I/O processes)
    """
    if min_io_processes > max_io_processes:
        raise ValueError(
            "Minimum I/O processes cannot be greater than maximum I/O processes."
        )
    if min_io_processes <= 0 or max_io_processes <= 0:
        raise ValueError("I/O processes must be positive integers.")

    best_allocation = (
        1,
        total_threads,
    )  # Start with all threads assigned to 1 worker if no better found
    max_workers = 0

    for io_processes in range(min_io_processes, max_io_processes + 1):
        if (
            total_threads >= io_processes
        ):  # Ensure there are enough threads to allocate at least these many I/O processes
            workers = total_threads // io_processes
            if (
                workers > max_workers
            ):  # Find the configuration with the maximum number of CPU workers
                max_workers = workers
                best_allocation = (workers, io_processes)

    return best_allocation
