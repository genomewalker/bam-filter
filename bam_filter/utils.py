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

log = logging.getLogger("my_logger")
log.setLevel(logging.INFO)
timestr = time.strftime("%Y%m%d-%H%M%S")


def is_debug():
    return logging.getLogger("my_logger").getEffectiveLevel() == logging.DEBUG


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
    value = int(val[:-1])

    if is_integer(value) & (unit in units) & (value > 0):
        if var == "--scale":
            if unit == "K":
                val = value * 1000
            elif unit == "M":
                val = value * 1000000
            elif unit == "G":
                val = value * 1000000000
            return val
        else:
            return val
    else:
        parser.error(
            "argument %s: Invalid value %s. Has to be an integer larger than 0 with the following suffix K, M or G"
            % (var, val)
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


defaults = {
    "min_read_length": 30,
    "min_read_count": 10,
    "min_expected_breadth_ratio": 0.5,
    "min_norm_entropy": "auto",
    "min_norm_gini": None,
    "min_avg_read_ani": 90.0,
    "min_read_ani": 90.0,
    "min_breadth": 0,
    "min_coverage_evenness": 0,
    "min_coverage_mean": 0,
    "prefix": None,
    "sort_memory": "1G",
    "reference_lengths": None,
    "scale": 1e6,
    "chunk_size": None,
    "plot": False,
}

help_msg = {
    "bam": "BAM file containing aligned reads",
    "threads": "Number of threads to use",
    "prefix": "Prefix used for the output files",
    "min_read_length": "Minimum read length",
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
    "min_coverage_mean": "Minimum coverage mean",
    "transform_cov_evenness": "Include those references that fullfi all filtering criteria but the coverage evenness is 0",
    "sort_memory": "Set maximum memory per thread for sorting; suffix K/M/G recognized",
    "scale": "Scale taxonomic abundance by this factor; suffix K/M recognized",
    "read_length_freqs": "Save a JSON file with the read length frequencies mapped to each reference",
    "read_hits_count": "Save a TSV file with the read hits frequencies mapped to each reference",
    "only_stats": "Only produce statistics and skip filtering",
    "only_stats_filtered": "Only filter statistics and skip BAM filtering",
    "plot": "Plot genome coverage plots",
    "sort_by_name": "Sort by read names",
    "chunk_size": "Chunk size for parallel processing",
    "help": "Help message",
    "debug": "Print debug messages",
    "reference_lengths": "File with references lengths",
    "version": "Print program version",
}


def get_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="A simple tool to calculate metrics from a BAM file and filter with uneven coverage.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # add subparser for filtering options:
    filter_args = parser.add_argument_group("filtering arguments")
    misc_args = parser.add_argument_group("miscellaneous arguments")
    parser.add_argument(
        "bam",
        type=lambda x: is_valid_file(parser, x, "bam"),
        help=help_msg["bam"],
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=lambda x: int(
            check_values(x, minval=1, maxval=1000, parser=parser, var="--threads")
        ),
        dest="threads",
        default=1,
        help=help_msg["threads"],
    )
    misc_args.add_argument(
        "--reference-trim-length",
        type=lambda x: int(
            check_values(
                x, minval=0, maxval=10000, parser=parser, var="---reference-trim-length"
            )
        ),
        dest="trim_ends",
        default=0,
        help=help_msg["trim_ends"],
    )
    misc_args.add_argument(
        "--trim-min",
        type=lambda x: int(
            check_values(x, minval=0, maxval=100, parser=parser, var="--trim-min")
        ),
        dest="trim_min",
        default=10,
        help=help_msg["trim_min"],
    )
    misc_args.add_argument(
        "--trim-max",
        type=lambda x: int(
            check_values(x, minval=0, maxval=100, parser=parser, var="--trim-max")
        ),
        dest="trim_max",
        default=90,
        help=help_msg["trim_max"],
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        default=defaults["prefix"],
        dest="prefix",
        help=help_msg["prefix"],
    )
    filter_args.add_argument(
        "-A",
        "--min-read-ani",
        type=lambda x: float(
            check_values(x, minval=0, maxval=100, parser=parser, var="--min-read-ani")
        ),
        default=defaults["min_read_ani"],
        dest="min_read_ani",
        help=help_msg["min_read_ani"],
    )
    filter_args.add_argument(
        "-l",
        "--min-read-length",
        type=lambda x: int(
            check_values(
                x, minval=1, maxval=100000, parser=parser, var="--min-read-length"
            )
        ),
        default=defaults["min_read_length"],
        dest="min_read_length",
        help=help_msg["min_read_length"],
    )
    filter_args.add_argument(
        "-n",
        "--min-read-count",
        type=lambda x: int(
            check_values(
                x, minval=1, maxval=100000, parser=parser, var="--min-read-count"
            )
        ),
        default=defaults["min_read_count"],
        dest="min_read_count",
        help=help_msg["min_read_count"],
    )
    filter_args.add_argument(
        "-b",
        "--min-expected-breadth-ratio",
        type=lambda x: float(
            check_values(
                x, minval=0, maxval=1, parser=parser, var="--min-expected-breadth-ratio"
            )
        ),
        default=defaults["min_expected_breadth_ratio"],
        dest="min_expected_breadth_ratio",
        help=help_msg["min_expected_breadth_ratio"],
    )
    filter_args.add_argument(
        "-e",
        "--min-normalized-entropy",
        type=lambda x: check_values_auto(
            x, minval=0, maxval=1, parser=parser, var="--min-normalized-entropy"
        ),
        default=defaults["min_norm_entropy"],
        dest="min_norm_entropy",
        help=help_msg["min_norm_entropy"],
    )
    filter_args.add_argument(
        "-g",
        "--min-normalized-gini",
        type=lambda x: check_values_auto(
            x, minval=0, maxval=1, parser=parser, var="--min-normalized-gini"
        ),
        default=defaults["min_norm_gini"],
        dest="min_norm_gini",
        help=help_msg["min_norm_gini"],
    )
    filter_args.add_argument(
        "-B",
        "--min-breadth",
        type=lambda x: float(
            check_values(x, minval=0, maxval=1, parser=parser, var="--min-breadth")
        ),
        default=defaults["min_breadth"],
        dest="min_breadth",
        help=help_msg["min_breadth"],
    )
    filter_args.add_argument(
        "-a",
        "--min-avg-read-ani",
        type=lambda x: float(
            check_values(
                x, minval=0, maxval=100, parser=parser, var="--min-avg-read-ani"
            )
        ),
        default=defaults["min_avg_read_ani"],
        dest="min_avg_read_ani",
        help=help_msg["min_avg_read_ani"],
    )
    filter_args.add_argument(
        "-c",
        "--min-coverage-evenness",
        type=lambda x: float(
            check_values(
                x, minval=0, maxval=1, parser=parser, var="--min-coverage-evenness"
            )
        ),
        default=defaults["min_coverage_evenness"],
        dest="min_coverage_evenness",
        help=help_msg["min_coverage_evenness"],
    )
    filter_args.add_argument(
        "-C",
        "--min-coverage-mean",
        type=lambda x: float(
            check_values(
                x, minval=0, maxval=1000000, parser=parser, var="--min-coverage-mean"
            )
        ),
        default=defaults["min_coverage_mean"],
        dest="min_coverage_mean",
        help=help_msg["min_coverage_mean"],
    )
    filter_args.add_argument(
        "--include-low-detection",
        dest="transform_cov_evenness",
        action="store_true",
        help=help_msg["transform_cov_evenness"],
    )
    parser.add_argument(
        "-m",
        "--sort-memory",
        type=lambda x: check_suffix(x, parser=parser, var="--sort-memory"),
        default=defaults["sort_memory"],
        dest="sort_memory",
        help=help_msg["sort_memory"],
    )
    parser.add_argument(
        "-N",
        "--sort-by-name",
        dest="sort_by_name",
        action="store_true",
        help=help_msg["sort_by_name"],
    )
    parser.add_argument(
        "--scale",
        type=lambda x: check_suffix(x, parser=parser, var="--scale"),
        default=defaults["scale"],
        dest="scale",
        help=help_msg["scale"],
    )
    # reference_lengths
    parser.add_argument(
        "-r",
        "--reference-lengths",
        type=lambda x: is_valid_file(parser, x, "reference_lengths"),
        default=defaults["reference_lengths"],
        dest="reference_lengths",
        help=help_msg["reference_lengths"],
    )
    parser.add_argument(
        "--read-length-freqs",
        dest="read_length_freqs",
        action="store_true",
        help=help_msg["read_length_freqs"],
    )
    parser.add_argument(
        "--read-hits-count",
        dest="read_hits_count",
        action="store_true",
        help=help_msg["read_hits_count"],
    )
    parser.add_argument(
        "--only-stats",
        dest="only_stats",
        action="store_true",
        help=help_msg["only_stats"],
    )
    parser.add_argument(
        "--plot",
        dest="plot",
        action="store_true",
        help=help_msg["plot"],
    )
    parser.add_argument(
        "--only-stats-filtered",
        dest="only_stats_filtered",
        action="store_true",
        help=help_msg["only_stats_filtered"],
    )
    parser.add_argument(
        "--chunk-size",
        type=lambda x: int(
            check_values(x, minval=1, maxval=100000, parser=parser, var="--chunk-size")
        ),
        default=defaults["chunk_size"],
        dest="chunk_size",
        help=help_msg["chunk_size"],
    )
    parser.add_argument(
        "--debug", dest="debug", action="store_true", help=help_msg["debug"]
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + __version__,
        help=help_msg["version"],
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


def do_parallel(parms, lst, func, threads):
    if is_debug():
        dfs = list(map(partial(func, parms=parms), lst))
    else:
        p = Pool(threads, initializer=initializer, initargs=(parms,))
        c_size = calc_chunksize(threads, len(lst))
        dfs = list(
            tqdm.tqdm(
                p.imap_unordered(
                    partial(func, parms=parms),
                    lst,
                    chunksize=c_size,
                ),
                total=len(lst),
                leave=False,
                ncols=80,
                desc=f"Components processed",
            )
        )
        p.close()
        p.join()
    return concat_df(dfs)


def do_parallel_lst(parms, lst, func, threads):
    if is_debug():
        lst = list(map(partial(func, parms=parms), lst))
    else:
        p = Pool(threads, initializer=initializer, initargs=(parms,))
        c_size = calc_chunksize(threads, len(lst))
        lst = list(
            tqdm.tqdm(
                p.imap_unordered(
                    partial(func, parms=parms),
                    lst,
                    chunksize=c_size,
                ),
                total=len(lst),
                leave=False,
                ncols=80,
                desc=f"Components processed",
            )
        )

    return lst


def get_components_large(parms, components, func, threads):
    dfs = list(
        tqdm.tqdm(
            map(partial(func, parms=parms), components),
            total=len(components),
            leave=False,
            ncols=80,
            desc=f"Components processed",
        )
    )
    return concat_df(dfs)


def clean_up(keep, temp_dir):
    if keep:
        logging.info(f"Cleaning up temporary files")
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


def create_output_files(prefix, bam):
    if prefix is None:
        prefix = bam.replace(".bam", "")
    # create output files
    out_files = {
        "stats": f"{prefix}_stats.tsv.gz",
        "stats_filtered": f"{prefix}_stats-filtered.tsv.gz",
        "bam_filtered_tmp": f"{prefix}.filtered.tmp.bam",
        "bam_filtered": f"{prefix}.filtered.bam",
        "read_length_freqs": f"{prefix}_read-length-freqs.json",
        "read_hits_count": f"{prefix}_read-hits-count.tsv.gz",
        "knee_plot": f"{prefix}_knee-plot.png",
        "coverage_plot_dir": f"{prefix}_coverage-plots",
    }
    return out_files
