import logging
import datatable as dt
import pandas as pd

# Import BAM utilities
from bam_filter.bam_utils import (
    check_bam_file,
    process_bam,
    filter_reference_BAM,
)

from bam_filter.utils import (
    get_arguments,
    create_output_files,
    concat_df,
    check_tmp_dir_exists,
    handle_warning,
    create_empty_output_files,
)
from bam_filter.entropy import find_knee
import json
import warnings
from collections import Counter
from functools import reduce
import os
import sys
from typing import Dict, List, TextIO, Union, Iterator, Any
from itertools import islice

log = logging.getLogger("my_logger")


class StreamingJSONWriter:
    """Helper class to write formatted JSON array with proper indentation."""

    def __init__(self, file: TextIO, indent: int = 4):
        self.file = file
        self.indent = indent
        self.first_item = True
        self._wrote_open = False

    def write_open(self):
        """Write the opening bracket of the JSON array."""
        self.file.write("[\n")
        self._wrote_open = True

    def write_close(self):
        """Write the closing bracket of the JSON array."""
        if self._wrote_open:
            self.file.write("\n]")

    def write_item(self, item: Any, obj_dict: Any = None):
        """Write a single item to the JSON array with proper formatting."""
        if not self._wrote_open:
            self.write_open()

        if not self.first_item:
            self.file.write(",\n")
        self.first_item = False

        json_str = json.dumps(item, default=obj_dict, ensure_ascii=False)
        indent_str = " " * self.indent
        self.file.write(f"{indent_str}{json_str}")


def write_lengths_to_file(
    data: List[Any], output_path: str, obj_dict: Any, batch_size: int = 10000
) -> None:
    """
    Write read lengths to file in a memory-efficient way with proper JSON formatting.

    Args:
        data: Input data containing items where second element is length
        output_path: Path to output file
        obj_dict: Function to handle object serialization
        batch_size: Size of batches for processing
    """
    logging.info("Writing read lengths to file...")

    with open(output_path, "w", encoding="utf-8") as outfile:
        writer = StreamingJSONWriter(outfile)
        data_iter = iter(data)

        while True:
            # Process data in batches
            batch = list(islice(data_iter, batch_size))
            if not batch:
                break

            # Extract valid lengths from current batch
            batch_lengths = [x[1] for x in batch if x[1] is not None]

            # Write each length
            for length in batch_lengths:
                writer.write_item(length, obj_dict)

        # Close the JSON array
        writer.write_close()

    logging.info(f"Wrote length list to {output_path}")


def obj_dict(obj):
    return obj.__dict__


def get_summary(obj):
    return obj.to_summary()


def get_lens(obj):
    return obj.get_read_length_freqs()


def filter_references(args):
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(levelname)s ::: %(asctime)s ::: %(message)s",
        datefmt="%H:%M:%S",
    )

    args = get_arguments()

    if args.debug:
        log.warning("Debug mode enabled")

    # Set datatable threads
    dt.options.nthreads = args.threads
    dt.options.progress.enabled = True
    dt.options.progress.clear_on_success = True

    tmp_dir = check_tmp_dir_exists(args.tmp_dir)
    log.info("Temporary directory: %s", tmp_dir.name)
    if args.trim_min >= args.trim_max:
        log.error("trim_min must be less than trim_max")
        exit(1)

    logging.getLogger("my_logger").setLevel(
        logging.DEBUG if args.debug else logging.INFO
    )

    if args.debug:
        warnings.showwarning = handle_warning
    else:
        warnings.filterwarnings("ignore")

    plot_coverage = False
    if args.coverage_plots is not None:
        plot_coverage = True

    if args.stats_filtered is None and (args.bam_filtered is not None):
        logging.error(
            "You need to specify a filtereds stats file to obtain the filtered BAM file"
        )
        exit(1)

    out_files = create_output_files(
        prefix=args.prefix,
        bam=args.bam,
        stats=args.stats,
        stats_filtered=args.stats_filtered,
        bam_filtered=args.bam_filtered,
        read_length_freqs=args.read_length_freqs,
        read_hits_count=args.read_hits_count,
        knee_plot=args.knee_plot,
        coverage_plots=args.coverage_plots,
        tmp_dir=tmp_dir,
        mode="filter",
    )

    bam = check_bam_file(
        bam=args.bam,
        threads=4,
        reference_lengths=args.reference_lengths,
        sort_memory=args.sort_memory,
        sorted_bam=out_files["sorted_bam"],
    )
    if bam is None:
        logging.warning("No reference sequences with alignments found in the BAM file")
        create_empty_output_files(out_files)
        sys.exit(0)

    data = process_bam(
        bam=bam,
        threads=args.threads,
        reference_lengths=args.reference_lengths,
        min_read_count=args.min_read_count,
        min_read_ani=args.min_read_ani,
        trim_ends=args.trim_ends,
        trim_min=args.trim_min,
        trim_max=args.trim_max,
        scale=args.scale,
        plot=plot_coverage,
        plots_dir=out_files["coverage_plot_dir"],
        read_length_freqs=args.read_length_freqs,
        output_files=out_files,
        sort_memory=args.sort_memory,
        low_memory=args.low_memory,
    )
    if args.low_memory:
        bam = out_files["bam_tmp_sorted"]

    logging.info("Reducing results to a single dataframe")
    data_df = [x[0] for x in data if x[0] is not None]
    data_df = [x for x in data_df if not x.empty]
    data_df = concat_df(data_df)

    if args.read_length_freqs is not None:
        write_lengths_to_file(
            data,
            out_files["read_length_freqs"],
            obj_dict,
            batch_size=1_000_000,  # Adjust based on available memory
        )

    if args.read_hits_count is not None:

        logging.info("Calculating read hits counts...")
        hits = [x[2] for x in data if x[2] is not None]

        # merge dicts and sum values
        hits = reduce(lambda x, y: x.update(y) or x, (Counter(dict(x)) for x in hits))

        # convert dict to datatable
        hits_dt = dt.Frame({"read": list(hits.keys()), "count": list(hits.values())})

        # Save using datatable
        hits_dt[:, :, dt.sort(-dt.f.count)].to_csv(
            out_files["read_hits_count"], sep="\t", header=True
        )

    logging.info(f"Writing reference statistics to {out_files['stats']}")
    # Convert pandas DataFrame to datatable
    data_dt = dt.Frame(data_df)
    data_dt.to_csv(out_files["stats"], sep="\t", header=True)

    if args.min_norm_entropy is None or args.min_norm_gini is None:
        filter_conditions = {
            "min_read_length": args.min_read_length,
            "max_read_length": args.max_read_length,
            "min_read_count": args.min_read_count,
            "min_expected_breadth_ratio": args.min_expected_breadth_ratio,
            "min_breadth": args.min_breadth,
            "min_avg_read_ani": args.min_avg_read_ani,
            "min_coverage_evenness": args.min_coverage_evenness,
            "min_coeff_var": args.min_coeff_var,
            "min_coverage_mean": args.min_coverage_mean,
        }
    elif args.min_norm_entropy == "auto" or args.min_norm_gini == "auto":
        if data_df.shape[0] > 1:
            min_norm_gini, min_norm_entropy = find_knee(
                data_df, out_plot_name=out_files["knee_plot"]
            )

            if min_norm_gini is None or min_norm_entropy is None:
                logging.warning(
                    "Could not find knee in entropy plot. Disabling filtering by entropy/gini inequality."
                )
                filter_conditions = {
                    "min_read_length": args.min_read_length,
                    "max_read_length": args.max_read_length,
                    "min_read_count": args.min_read_count,
                    "min_expected_breadth_ratio": args.min_expected_breadth_ratio,
                    "min_breadth": args.min_breadth,
                    "min_avg_read_ani": args.min_avg_read_ani,
                    "min_coverage_evenness": args.min_coverage_evenness,
                    "min_coeff_var": args.min_coeff_var,
                    "min_coverage_mean": args.min_coverage_mean,
                }
            else:
                filter_conditions = {
                    "min_read_length": args.min_read_length,
                    "max_read_length": args.max_read_length,
                    "min_read_count": args.min_read_count,
                    "min_expected_breadth_ratio": args.min_expected_breadth_ratio,
                    "min_breadth": args.min_breadth,
                    "min_avg_read_ani": args.min_avg_read_ani,
                    "min_coverage_evenness": args.min_coverage_evenness,
                    "min_coeff_var": args.min_coeff_var,
                    "min_coverage_mean": args.min_coverage_mean,
                    "min_norm_entropy": min_norm_entropy,
                    "min_norm_gini": min_norm_gini,
                }
        else:
            min_norm_gini = 0.4
            min_norm_entropy = 0.75
            logging.warning(
                f"There's only one genome. Using min_norm_gini={min_norm_gini} and min_norm_entropy={min_norm_entropy}. Please check the results."
            )
            filter_conditions = {
                "min_read_length": args.min_read_length,
                "max_read_length": args.max_read_length,
                "min_read_count": args.min_read_count,
                "min_expected_breadth_ratio": args.min_expected_breadth_ratio,
                "min_breadth": args.min_breadth,
                "min_avg_read_ani": args.min_avg_read_ani,
                "min_coverage_evenness": args.min_coverage_evenness,
                "min_coeff_var": args.min_coeff_var,
                "min_coverage_mean": args.min_coverage_mean,
                "min_norm_entropy": min_norm_entropy,
                "min_norm_gini": min_norm_gini,
            }
    else:
        min_norm_gini, min_norm_entropy = args.min_norm_gini, args.min_norm_entropy
        filter_conditions = {
            "min_read_length": args.min_read_length,
            "max_read_length": args.max_read_length,
            "min_read_count": args.min_read_count,
            "min_expected_breadth_ratio": args.min_expected_breadth_ratio,
            "min_breadth": args.min_breadth,
            "min_avg_read_ani": args.min_avg_read_ani,
            "min_coverage_evenness": args.min_coverage_evenness,
            "min_coeff_var": args.min_coeff_var,
            "min_coverage_mean": args.min_coverage_mean,
            "min_norm_entropy": min_norm_entropy,
            "min_norm_gini": min_norm_gini,
        }

    if args.stats_filtered is not None:
        filter_reference_BAM(
            bam=bam,
            df=data_df,
            filter_conditions=filter_conditions,
            transform_cov_evenness=args.transform_cov_evenness,
            threads=args.threads,
            out_files=out_files,
            sort_memory=args.sort_memory,
            sort_by_name=args.sort_by_name,
            min_read_ani=args.min_read_ani,
            disable_sort=args.disable_sort,
            write_bam=args.bam_filtered is not None,
        )
    else:
        logging.info("Skipping filtering of reference BAM file.")

    if args.low_memory:
        os.remove(out_files["bam_tmp_sorted"])
    logging.info("ALL DONE.")
