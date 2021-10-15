import argparse
import sys
import gzip
import os
import pathlib
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
import pathlib
import uuid
import subprocess
from itertools import chain
from statistics import mean

log = logging.getLogger("my_logger")
log.setLevel(logging.INFO)
timestr = time.strftime("%Y%m%d-%H%M%S")


def is_debug():
    return logging.getLogger("my_logger").getEffectiveLevel() == logging.DEBUG


def check_values(val, minval, maxval, parser, var):
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
    return value


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
        if os.path.isfile(arg):
            parser.error("argument %s: The file %s does not exist!" % (var, arg))
        else:
            parser.error("argument %s: The directory %s does not exist!" % (var, arg))
    else:
        if os.path.isfile(arg):
            return get_open_func(arg)(arg, "rt")  # return an open file handle
        else:
            return arg


defaults = {}

help_msg = {
    "bam": "BAM file containing aligned reads",
    "threads": "Number of threads to use",
    "output": "Output file name",
    "breadthExpRatio": "Breadth to expected breadth ratio",
    "covEvenness": "Eveness of coverage",
    "nReads": "Number of reads",
    "aniNM": "ANI between NM and MD",
    "aniMD": "ANI between MD and NM",
    "help": "Help message",
}


def get_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="Calculate metrics for a BAM file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
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
        default=1,
        help=help_msg["threads"],
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="metrics.tsv",
        help=help_msg["output"],
    )
    parser.add_argument(
        "-B",
        "--breadthExpRatio",
        type=lambda x: check_values(
            x, minval=0, maxval=1, parser=parser, var="--breadthExpRatio"
        ),
        default=0.0,
        help=help_msg["breadthExpRatio"],
    )
    parser.add_argument(
        "-C",
        "--covEvenness",
        type=lambda x: check_values(
            x, minval=0, maxval=1, parser=parser, var="--covEvenness"
        ),
        default=0.0,
        help=help_msg["covEvenness"],
    )
    # nReads
    parser.add_argument(
        "-n",
        "--nReads",
        type=lambda x: check_values(
            x, minval=0, maxval=1e9, parser=parser, var="--nReads"
        ),
        default=0.0,
        help=help_msg["nReads"],
    )
    # aniNM
    parser.add_argument(
        "-a",
        "--aniNM",
        type=lambda x: check_values(
            x, minval=0, maxval=100, parser=parser, var="--aniNM"
        ),
        default=0.0,
        help=help_msg["aniNM"],
    )
    # aniMD
    parser.add_argument(
        "-A",
        "--aniMD",
        type=lambda x: check_values(
            x, minval=0, maxval=100, parser=parser, var="--aniMD"
        ),
        default=0.0,
        help=help_msg["aniMD"],
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
