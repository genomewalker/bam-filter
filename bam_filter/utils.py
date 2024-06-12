import argparse
import sys
import gzip
import os
import shutil
import logging
import pandas as pd
from multiprocessing import Pool
from functools import partial
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from os import devnull
import tqdm
from bam_filter import __version__
import time
from itertools import chain
import numpy as np
from pathlib import Path
import pysam
import tempfile

log = logging.getLogger("my_logger")
log.setLevel(logging.INFO)
timestr = time.strftime("%Y%m%d-%H%M%S")


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
            log.error(f"Temporary directory {tmpdir} does not exist")
            exit(1)
        tmpdir = tempfile.TemporaryDirectory(dir=os.path.abspath(tmpdir))
        # Check if tmpdir has more than 107 characters
        if len(tmpdir.name) > 107:
            log.error(f"Temporary directory {tmpdir.name} has more than 107 characters")
            exit(1)
    return tmpdir


def is_debug():
    return logging.getLogger("my_logger").getEffectiveLevel() == logging.DEBUG


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
    if mode == "weight":
        target_weight = scale * max(input_dict.values())
        total_weight = sum(input_dict.values())
        target_weight = max(
            scale * max(input_dict.values()), max_entries_per_chunk or 0
        )
        num_chunks = max(num_cores, (total_weight // target_weight) + 1)
    elif mode == "entries":
        # split that it has at least 1000 entries per chunk
        total_entries = len(input_dict)
        entries = list(input_dict.keys())
        num_chunks = max(num_cores, (total_entries // num_entries) + 1)
        num_chunks = 100
        # split entries into num_chunks
        chunks = [entries[i::num_chunks] for i in range(num_chunks)]

    # Refine chunks
    if mode == "weight":
        sorted_keys = sorted(input_dict, key=input_dict.get, reverse=True)
        chunks = [[] for _ in range(num_chunks)]
        chunk_weights = [0] * num_chunks

        for key in tqdm.tqdm(
            sorted_keys,
            desc="Distributing keys",
            unit_scale=True,
            unit_divisor=500,
            leave=False,
            ncols=80,
            unit=" keys",
        ):
            min_chunk_idx = chunk_weights.index(min(chunk_weights))
            chunks[min_chunk_idx].append(key)
            chunk_weights[min_chunk_idx] += input_dict[key]

        for _ in range(refinement_steps):
            initial_balance = max(len(chunk) for chunk in chunks) - min(
                len(chunk) for chunk in chunks
            )
            chunks = refine_chunks(chunks, input_dict, chunk_weights, target_weight)
            current_balance = max(len(chunk) for chunk in chunks) - min(
                len(chunk) for chunk in chunks
            )
            if current_balance >= initial_balance:
                break  # No improvement, exit the loop

    if verbose:
        print_chunk_stats(chunks, input_dict)

    return chunks


def print_chunk_stats(chunks, input_dict):
    for i, chunk in enumerate(chunks, 1):
        weights = [input_dict[key] for key in chunk]
        print(
            f"Chunk {i}: Total = {sum(weights)}, Min = {min(weights)}, Max = {max(weights)}, Avg = {sum(weights) / len(chunk)}"
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


# From: https://note.nkmk.me/en/python-check-int-float/
def is_integer(n):
    try:
        float(n)
    except ValueError:
        return False
    else:
        return float(n).is_integer()


# function to check if the input value has K, M or G suffix in it


def check_suffix(val, parser, var):
    if var == "--scale":
        units = ["K", "M"]
    else:
        units = ["K", "M", "G"]

    unit = val[-1]
    value = val[:-1]

    if is_integer(value) and (unit in units) and (int(value) > 0):
        value = int(value)
        if var == "--scale":
            if unit == "K":
                val = value * 1000
            elif unit == "M":
                val = value * 1000000
            return str(val)
        else:
            if unit == "K":
                val = value * 1024
            elif unit == "M":
                val = value * 1024 * 1024
            elif unit == "G":
                val = value * 1024 * 1024 * 1024
            return val
    else:
        parser.error(
            "argument %s: Invalid value %s. Has to be an integer larger than 0 with the following suffix K, M or G"
            % (var, val)
        )


# Example usage with argparse
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process some integers.")
    parser.add_argument(
        "--scale",
        type=lambda x: check_suffix(x, parser, "--scale"),
        help="Scale value with K or M suffix",
    )
    parser.add_argument(
        "--memory",
        type=lambda x: check_suffix(x, parser, "--memory"),
        help="Memory value with K, M, or G suffix",
    )

    args = parser.parse_args()
    print("Scale:", args.scale)
    print("Memory:", args.memory)


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
    "max_read_length": np.Inf,
    "min_read_count": 3,
    "min_expected_breadth_ratio": 0,
    "min_norm_entropy": 0,
    "min_norm_gini": 1.0,
    "min_avg_read_ani": 90.0,
    "min_read_ani": 90.0,
    "min_breadth": 0,
    "min_coverage_evenness": 0,
    "min_coeff_var": np.Inf,
    "min_coverage_mean": 0,
    "prefix": None,
    "sort_memory": "1G",
    "reference_lengths": None,
    "scale": 1e6,
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
    "reassign_iters": 25,
    "reassign_scale": 0.9,
    "reassign_match_reward": 1,
    "reassign_mismatch_penalty": -2,
    "reassign_gap_open_penalty": 5,
    "reassign_gap_extension_penalty": 2,
    "reassign_lambda": 1.33,
    "reassign_k": 0.621,
    "rank_lca": "species",
    "lca_summary": None,
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
    "min_norm_entropy": "Minimum normalized entropy",
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
    "debug": "Print debug messages",
    "reference_lengths": "File with references lengths",
    "low_memory": "Activate the low memory mode",
    "reassign": "Run an EM algorithm to reassign reads to references",
    "reassign_method": "Method for the EM algorithm",
    "reassign_iters": "Number of iterations for the EM algorithm",
    "reassign_scale": "Scale to select the best weithing alignments",
    "reassign_match_reward": "Match reward for the alignment score ",
    "reassign_mismatch_penalty": "Mismatch penalty for the alignment score ",
    "reassign_gap_open_penalty": "Gap open penalty for alignment score computation",
    "reassign_gap_extension_penalty": "Gap extension penalty for the alignment score",
    "reassign_lambda": "Lambda parameter for the alignment score",
    "reassign_k": "K parameter for the alignment score algorithm",
    "lca": "Calculate LCA for each read and estimate abundances",
    "names": "Names dmp file from taxonomy",
    "nodes": "Nodes dmp file from taxonomy",
    "acc2taxid": "acc2taxid file from taxonomy",
    "rank_lca": "Rank to use for LCA calculation",
    "lca_summary": "Save a TSV file with the LCA summary",
    "lca_missing": "Save a TSV file with references with missing taxids",
    "lca_stats": "A TSV file from the filter subcommand",
    "custom": "Use custom taxdump files",
    "version": "Print program version",
    "max_memory": "Maximum memory to use for the EM algorithm",
}


def get_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="A simple tool to calculate metrics from a BAM file and filter with uneven coverage.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + __version__,
        help=help_msg["version"],
    )
    parser.add_argument(
        "--debug", dest="debug", action="store_true", help=help_msg["debug"]
    )

    sub_parsers = parser.add_subparsers(
        help="positional arguments",
        dest="action",
    )

    # Create parent subparser. Note `add_help=False` and creation via `argparse.`
    parent_parser = argparse.ArgumentParser(add_help=False)

    required = parent_parser.add_argument_group("required arguments")
    required.add_argument(
        "--bam",
        required=True,
        dest="bam",
        type=lambda x: is_valid_file(parser, x, "bam"),
        help=help_msg["bam"],
    )
    # required = parent_parser.add_argument_group("required arguments")
    optional = parent_parser.add_argument_group("optional arguments")
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
        dest="reference_lengths",
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
    # Create the parser sub-command for db creation
    parser_reassign = sub_parsers.add_parser(
        "reassign",
        help="Reassign reads to references using an EM algorithm",
        parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # create the parser sub-commands
    parser_filter = sub_parsers.add_parser(
        "filter",
        help="Filter references based on coverage and other metrics",
        parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_lca = sub_parsers.add_parser(
        "lca",
        help="Calculate LCA for each read and estimate abundances at each rank",
        parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # createdb_required_args = parser_createdb.add_argument_group("required arguments")
    # reassign_required_args = parser_reassign.add_argument_group(
    #     "Re-assign required arguments"
    # )
    reassign_optional_args = parser_reassign.add_argument_group(
        "Re-assign optional arguments"
    )
    misc_reassign_args = parser_reassign.add_argument_group("miscellaneous arguments")

    filter_required_args = parser_filter.add_argument_group("Filter required arguments")
    # filter_optional_args = parser_filter.add_argument_group("Filter optional arguments")

    # lca_required_args = parser_lca.add_argument_group("LCA required arguments")
    lca_optional_args = parser_lca.add_argument_group("LCA optional arguments")

    # add subparser for filtering options:
    # reassign_args = parser.add_argument_group("reassign arguments")
    filtering_filt_args = parser_filter.add_argument_group("filtering arguments")
    # lca_args = parser.add_argument_group("lca arguments")
    misc_filter_args = parser_filter.add_argument_group("miscellaneous arguments")
    out_filter_args = parser_filter.add_argument_group("output arguments")
    # parser.add_argument(
    #     "--bam",
    #     required=True,
    #     dest="bam",
    #     type=lambda x: is_valid_file(parser, x, "bam"),
    #     help=help_msg["bam"],
    # )
    # parser.add_argument(
    #     "-t",
    #     "--threads",
    #     type=lambda x: int(
    #         check_values(x, minval=1, maxval=1000, parser=parser, var="--threads")
    #     ),
    #     dest="threads",
    #     metavar="INT",
    #     default=1,
    #     help=help_msg["threads"],
    # )

    reassign_optional_args.add_argument(
        "-i",
        "--iters",
        type=lambda x: int(
            check_values(
                x, minval=0, maxval=100000, parser=parser, var="--reassign-n-iters"
            )
        ),
        metavar="INT",
        default=defaults["reassign_iters"],
        dest="reassign_iters",
        help=help_msg["reassign_iters"],
    )
    reassign_optional_args.add_argument(
        "-s",
        "--scale",
        type=lambda x: float(
            check_values(x, minval=0, maxval=1, parser=parser, var="--scale")
        ),
        metavar="FLOAT",
        default=defaults["reassign_scale"],
        dest="reassign_scale",
        help=help_msg["reassign_scale"],
    )
    reassign_optional_args.add_argument(
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
    reassign_optional_args.add_argument(
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
    reassign_optional_args.add_argument(
        "-L",
        "--max-read-length",
        type=lambda x: int(
            check_values(
                x, minval=1, maxval=np.Inf, parser=parser, var="--max-read-length"
            )
        ),
        default=defaults["max_read_length"],
        metavar="INT",
        dest="max_read_length",
        help=help_msg["max_read_length"],
    )
    reassign_optional_args.add_argument(
        "-n",
        "--min-read-count",
        type=lambda x: int(
            check_values(
                x, minval=1, maxval=np.Inf, parser=parser, var="--min-read-count"
            )
        ),
        default=defaults["min_read_count"],
        metavar="INT",
        dest="min_read_count",
        help=help_msg["min_read_count"],
    )
    reassign_optional_args.add_argument(
        "--match-reward",
        type=lambda x: int(
            check_values(
                x, minval=0, maxval=np.Inf, parser=parser, var="--match-reward"
            )
        ),
        default=defaults["reassign_match_reward"],
        metavar="INT",
        dest="match_reward",
        help=help_msg["reassign_match_reward"],
    )
    reassign_optional_args.add_argument(
        "--mismatch-penalty",
        type=lambda x: int(
            check_values(
                x, minval=-np.Inf, maxval=0, parser=parser, var="--mismatch-penalty"
            )
        ),
        default=defaults["reassign_mismatch_penalty"],
        metavar="INT",
        dest="mismatch_penalty",
        help=help_msg["reassign_mismatch_penalty"],
    )
    reassign_optional_args.add_argument(
        "--gap-open-penalty",
        type=lambda x: int(
            check_values(
                x, minval=0, maxval=np.Inf, parser=parser, var="--gap-open-penalty"
            )
        ),
        default=defaults["reassign_gap_open_penalty"],
        metavar="INT",
        dest="gap_open_penalty",
        help=help_msg["reassign_gap_open_penalty"],
    )
    reassign_optional_args.add_argument(
        "--gap-extension-penalty",
        type=lambda x: int(
            check_values(
                x, minval=0, maxval=np.Inf, parser=parser, var="--gap-extension-penalty"
            )
        ),
        default=defaults["reassign_gap_extension_penalty"],
        metavar="INT",
        dest="gap_extension_penalty",
        help=help_msg["reassign_gap_extension_penalty"],
    )
    reassign_optional_args.add_argument(
        "--lambda",
        type=lambda x: float(
            check_values(x, minval=0, maxval=np.Inf, parser=parser, var="--lambda")
        ),
        default=defaults["reassign_lambda"],
        metavar="FLOAT",
        dest="lambda_value",
        help=help_msg["reassign_lambda"],
    )
    reassign_optional_args.add_argument(
        "-k",
        type=lambda x: float(
            check_values(x, minval=0, maxval=np.Inf, parser=parser, var="-K")
        ),
        default=defaults["reassign_k"],
        metavar="FLOAT",
        dest="K_value",
        help=help_msg["reassign_k"],
    )
    reassign_optional_args.add_argument(
        "-o",
        "--out-bam",
        dest="bam_reassigned",
        default=defaults["bam_reassigned"],
        metavar="FILE",
        type=str,
        nargs="?",
        const="",
        help=help_msg["bam_reassigned"],
    )
    reassign_optional_args.add_argument(
        "-m",
        "--sort-memory",
        type=lambda x: check_suffix(x, parser=parser, var="--sort-memory"),
        default=defaults["sort_memory"],
        metavar="STR",
        dest="sort_memory",
        help=help_msg["sort_memory"],
    )
    reassign_optional_args.add_argument(
        "-M",
        "--max-memory",
        type=lambda x: check_suffix(x, parser=parser, var="--max-memory"),
        default=None,
        metavar="INT",
        dest="max_memory",
        help=help_msg["max_memory"],
    )
    reassign_optional_args.add_argument(
        "-N",
        "--sort-by-name",
        dest="sort_by_name",
        action="store_true",
        help=help_msg["sort_by_name"],
    )
    reassign_optional_args.add_argument(
        "--tmp-dir",
        type=str,
        default=defaults["tmp_dir"],
        metavar="DIR",
        dest="tmp_dir",
        help=help_msg["tmp_dir"],
    )
    misc_reassign_args.add_argument(
        "--disable-sort",
        dest="disable_sort",
        action="store_true",
        help=help_msg["disable_sort"],
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
    filtering_filt_args.add_argument(
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
        help=help_msg["min_read_length"],
    )
    filtering_filt_args.add_argument(
        "-L",
        "--max-read-length",
        type=lambda x: int(
            check_values(
                x, minval=1, maxval=np.Inf, parser=parser, var="--max-read-length"
            )
        ),
        default=defaults["max_read_length"],
        metavar="INT",
        dest="max_read_length",
        help=help_msg["max_read_length"],
    )
    filtering_filt_args.add_argument(
        "-n",
        "--min-read-count",
        type=lambda x: int(
            check_values(
                x, minval=1, maxval=np.Inf, parser=parser, var="--min-read-count"
            )
        ),
        default=defaults["min_read_count"],
        metavar="INT",
        dest="min_read_count",
        help=help_msg["min_read_count"],
    )
    filtering_filt_args.add_argument(
        "-b",
        "--min-expected-breadth-ratio",
        type=lambda x: float(
            check_values(
                x, minval=0, maxval=1, parser=parser, var="--min-expected-breadth-ratio"
            )
        ),
        metavar="FLOAT",
        default=defaults["min_expected_breadth_ratio"],
        dest="min_expected_breadth_ratio",
        help=help_msg["min_expected_breadth_ratio"],
    )
    filtering_filt_args.add_argument(
        "-e",
        "--min-normalized-entropy",
        type=lambda x: check_values_auto(
            x, minval=0, maxval=1, parser=parser, var="--min-normalized-entropy"
        ),
        default=defaults["min_norm_entropy"],
        metavar="FLOAT",
        dest="min_norm_entropy",
        help=help_msg["min_norm_entropy"],
    )
    filtering_filt_args.add_argument(
        "-g",
        "--min-normalized-gini",
        type=lambda x: check_values_auto(
            x, minval=0, maxval=1, parser=parser, var="--min-normalized-gini"
        ),
        default=defaults["min_norm_gini"],
        metavar="FLOAT",
        dest="min_norm_gini",
        help=help_msg["min_norm_gini"],
    )
    filtering_filt_args.add_argument(
        "-B",
        "--min-breadth",
        type=lambda x: float(
            check_values(x, minval=0, maxval=1, parser=parser, var="--min-breadth")
        ),
        default=defaults["min_breadth"],
        metavar="FLOAT",
        dest="min_breadth",
        help=help_msg["min_breadth"],
    )
    filtering_filt_args.add_argument(
        "-a",
        "--min-avg-read-ani",
        type=lambda x: float(
            check_values(
                x, minval=0, maxval=100, parser=parser, var="--min-avg-read-ani"
            )
        ),
        metavar="FLOAT",
        default=defaults["min_avg_read_ani"],
        dest="min_avg_read_ani",
        help=help_msg["min_avg_read_ani"],
    )
    filtering_filt_args.add_argument(
        "-c",
        "--min-coverage-evenness",
        type=lambda x: float(
            check_values(
                x, minval=0, maxval=1, parser=parser, var="--min-coverage-evenness"
            )
        ),
        metavar="FLOAT",
        default=defaults["min_coverage_evenness"],
        dest="min_coverage_evenness",
        help=help_msg["min_coverage_evenness"],
    )
    filtering_filt_args.add_argument(
        "-V",
        "--min-coeff-var",
        type=lambda x: float(
            check_values(
                x, minval=0, maxval=np.Inf, parser=parser, var="--min-evenness"
            )
        ),
        default=defaults["min_coeff_var"],
        metavar="FLOAT",
        dest="min_coeff_var",
        help=help_msg["min_coeff_var"],
    )
    filtering_filt_args.add_argument(
        "-C",
        "--min-coverage-mean",
        type=lambda x: float(
            check_values(
                x, minval=0, maxval=1000000, parser=parser, var="--min-coverage-mean"
            )
        ),
        default=defaults["min_coverage_mean"],
        metavar="FLOAT",
        dest="min_coverage_mean",
        help=help_msg["min_coverage_mean"],
    )
    filtering_filt_args.add_argument(
        "--include-low-detection",
        dest="transform_cov_evenness",
        action="store_true",
        help=help_msg["transform_cov_evenness"],
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
    # reference_lengths
    # filter_optional_args.add_argument(
    #     "-r",
    #     "--reference-lengths",
    #     type=lambda x: is_valid_file(parser, x, "reference_lengths"),
    #     metavar="FILE",
    #     default=defaults["reference_lengths"],
    #     dest="reference_lengths",
    #     help=help_msg["reference_lengths"],
    # )
    filter_required_args.add_argument(
        "--stats",
        dest="stats",
        default=defaults["stats"],
        type=str,
        metavar="FILE",
        nargs="?",
        const="",
        required=True,
        help=help_msg["stats"],
    )
    out_filter_args.add_argument(
        "--stats-filtered",
        dest="stats_filtered",
        default=defaults["stats_filtered"],
        type=str,
        metavar="FILE",
        nargs="?",
        const="",
        help=help_msg["stats_filtered"],
    )
    out_filter_args.add_argument(
        "--bam-filtered",
        dest="bam_filtered",
        default=defaults["bam_filtered"],
        metavar="FILE",
        type=str,
        nargs="?",
        const="",
        help=help_msg["bam_filtered"],
    )
    out_filter_args.add_argument(
        "--read-length-freqs",
        dest="read_length_freqs",
        default=defaults["read_length_freqs"],
        metavar="FILE",
        type=str,
        nargs="?",
        const="",
        help=help_msg["read_length_freqs"],
    )
    out_filter_args.add_argument(
        "--read-hits-count",
        dest="read_hits_count",
        default=defaults["read_hits_count"],
        metavar="FILE",
        type=str,
        nargs="?",
        const="",
        help=help_msg["read_hits_count"],
    )
    out_filter_args.add_argument(
        "--knee-plot",
        dest="knee_plot",
        default=defaults["knee_plot"],
        metavar="FILE",
        type=str,
        nargs="?",
        const="",
        help=help_msg["knee_plot"],
    )
    out_filter_args.add_argument(
        "--coverage-plots",
        dest="coverage_plots",
        metavar="FILE",
        default=defaults["coverage_plots"],
        type=str,
        nargs="?",
        const="",
        help=help_msg["coverage_plots"],
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
    misc_filter_args.add_argument(
        "--tmp-dir",
        type=str,
        default=defaults["tmp_dir"],
        metavar="DIR",
        dest="tmp_dir",
        help=help_msg["tmp_dir"],
    )
    misc_filter_args.add_argument(
        "--low-memory",
        dest="low_memory",
        action="store_true",
        help=help_msg["low_memory"],
    )

    lca_optional_args.add_argument(
        "--names",
        metavar="FILE",
        type=lambda x: is_valid_file(parser, x, "names"),
        dest="names",
        help=help_msg["names"],
    )
    lca_optional_args.add_argument(
        "--nodes",
        metavar="FILE",
        type=lambda x: is_valid_file(parser, x, "nodes"),
        dest="nodes",
        help=help_msg["nodes"],
    )
    lca_optional_args.add_argument(
        "--acc2taxid",
        metavar="FILE",
        type=lambda x: is_valid_file(parser, x, "acc2taxid"),
        dest="acc2taxid",
        help=help_msg["acc2taxid"],
    )
    lca_optional_args.add_argument(
        "--lca-rank",
        metavar="STR",
        type=lambda x: str(check_lca_ranks(x, parser=parser, var="--lca-rank")),
        default=defaults["rank_lca"],
        dest="rank_lca",
        help=help_msg["rank_lca"],
    )
    lca_optional_args.add_argument(
        "--lca-summary",
        dest="lca_summary",
        metavar="FILE",
        default=defaults["lca_summary"],
        type=str,
        nargs="?",
        const="",
        help=help_msg["lca_summary"],
    )
    lca_optional_args.add_argument(
        "--scale",
        type=lambda x: check_suffix(x, parser=parser, var="--scale"),
        default=defaults["scale"],
        dest="scale",
        metavar="STR",
        help=help_msg["scale"],
    )
    lca_optional_args.add_argument(
        "-m",
        "--sort-memory",
        type=lambda x: check_suffix(x, parser=parser, var="--sort-memory"),
        default=defaults["sort_memory"],
        metavar="STR",
        dest="sort_memory",
        help=help_msg["sort_memory"],
    )
    lca_optional_args.add_argument(
        "--custom",
        dest="custom",
        action="store_true",
        help=help_msg["custom"],
    )
    lca_optional_args.add_argument(
        "--stats",
        dest="lca_stats",
        default=defaults["stats"],
        type=str,
        metavar="FILE",
        nargs="?",
        const="",
        required=False,
        help=help_msg["lca_stats"],
    )
    args = parser.parse_args(None if sys.argv[1:] else ["-h"])
    return args


@contextmanager
def suppress_stdout():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(devnull, "w") as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)


def applyParallel(dfGrouped, func, threads, parms):
    p = Pool(threads)
    func = partial(func, parms=parms)
    ret_list = tqdm.tqdm(
        p.map(func, [group for name, group in dfGrouped]),
        total=len([group for name, group in dfGrouped]),
    )
    p.close()
    p.join()
    return pd.concat(ret_list)


def fast_flatten(input_list):
    return list(chain.from_iterable(input_list))


def concat_df(frames):
    COLUMN_NAMES = frames[0].columns
    df_dict = dict.fromkeys(COLUMN_NAMES, [])
    for col in COLUMN_NAMES:
        extracted = (frame[col] for frame in frames)
        # Flatten and save to df_dict
        df_dict[col] = fast_flatten(extracted)
    df = pd.DataFrame.from_dict(df_dict)[COLUMN_NAMES]
    return df


def initializer(init_data):
    global parms
    parms = init_data


def clean_up(keep, temp_dir):
    if keep:
        logging.info("Cleaning up temporary files")
        logging.shutdown()
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

    if tmp_dir is not None:
        tmp_dir = tmp_dir.name
    else:
        tmp_dir = check_tmp_dir_exists(tmp_dir).name

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
            "bam_filtered_tmp": f"{tmp_dir}/{prefix}.filtered.tmp.bam",
            "bam_filtered": bam_filtered,
            "read_length_freqs": read_length_freqs,
            "read_hits_count": read_hits_count,
            "knee_plot": knee_plot,
            "coverage_plot_dir": coverage_plots,
            "bam_tmp": f"{tmp_dir}/{prefix}.tmp.bam",
            "bam_tmp_sorted": f"{tmp_dir}/{prefix}.tmp.sorted.bam",
        }
    elif mode == "reassign":
        out_files = {
            "bam_reassigned_tmp": f"{tmp_dir}/{prefix}.reassigned.tmp.bam",
            "bam_reassigned_sorted": f"{tmp_dir}/{prefix}.reassigned.sorted.bam",
            "bam_reassigned": bam_reassigned,
        }
    elif mode == "lca":
        out_files = {
            "lca_summary": lca_summary,
        }
    else:
        log.error("Mode not recognized")
        exit(1)
    out_files["tmp_dir"] = tmp_dir
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
