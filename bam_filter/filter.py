"""
This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not,
see <https://www.gnu.org/licenses/>.
"""

import logging
import pandas as pd

# from bam_filter.sam_utils import process_bam, filter_reference_BAM, check_bam_file
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

log = logging.getLogger("my_logger")


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
    # data = list(filter(None, data))
    data_df = [x[0] for x in data if x[0] is not None]
    # remove empty dataframes
    data_df = [x for x in data_df if not x.empty]
    data_df = concat_df(data_df)
    if args.read_length_freqs is not None:
        logging.info("Calculating read length frequencies...")
        lens = [x[1] for x in data if x[1] is not None]
        lens = json.dumps(lens, default=obj_dict, ensure_ascii=False, indent=4)
        with open(out_files["read_length_freqs"], "w", encoding="utf-8") as outfile:
            print(lens, file=outfile)

    if args.read_hits_count is not None:
        logging.info("Calculating read hits counts...")
        hits = [x[2] for x in data if x[2] is not None]

        # merge dicts and sum values
        hits = reduce(lambda x, y: x.update(y) or x, (Counter(dict(x)) for x in hits))
        # hits = sum(map(Counter, hits), Counter())

        # convert dict to dataframe
        hits = (
            pd.DataFrame.from_dict(hits, orient="index", columns=["count"])
            .rename_axis("read")
            .reset_index()
            .sort_values(by="count", ascending=False)
        )

        hits.to_csv(
            out_files["read_hits_count"],
            sep="\t",
            index=False,
        )

    logging.info(f"Writing reference statistics to {out_files['stats']}")
    data_df.to_csv(out_files["stats"], sep="\t", index=False)

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
            min_norm_gini = 0.5
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
    # check if sorted BAM file exists, if yes remove it
    logging.info("ALL DONE.")
