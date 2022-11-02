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
from bam_filter.sam_utils import process_bam, filter_reference_BAM
from bam_filter.utils import get_arguments, create_output_files
from bam_filter.entropy import find_knee
import json
import warnings

log = logging.getLogger("my_logger")


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


def obj_dict(obj):
    return obj.__dict__


def main():

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(levelname)s ::: %(asctime)s ::: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    args = get_arguments()
    logging.getLogger("my_logger").setLevel(
        logging.DEBUG if args.debug else logging.INFO
    )

    if args.debug:
        warnings.showwarning = handle_warning
    else:
        warnings.filterwarnings("ignore")

    out_files = create_output_files(prefix=args.prefix, bam=args.bam)

    data = process_bam(
        bam=args.bam,
        threads=args.threads,
        reference_lengths=args.reference_lengths,
        min_read_count=args.min_read_count,
        scale=args.scale,
        sort_memory=args.sort_memory,
        plot=args.plot,
        plots_dir=out_files["coverage_plot_dir"],
        chunksize=args.chunk_size,
    )

    data = list(filter(None, data))

    data_df = pd.DataFrame([x.to_summary() for x in data])

    if args.read_length_freqs:
        lens = [x.get_read_length_freqs() for x in data]
        lens = json.dumps(lens, default=obj_dict, ensure_ascii=False, indent=4)
        with open(out_files["read_length_freqs"], "w", encoding="utf-8") as outfile:
            print(lens, file=outfile)

    logging.info(f"Writing reference statistics to {out_files['stats']}")
    data_df.to_csv(out_files["stats"], sep="\t", index=False, compression="gzip")

    if args.min_norm_entropy == "auto" or args.min_norm_gini == "auto":
        if data_df.shape[0] > 1:
            min_norm_gini, min_norm_entropy = find_knee(
                data_df, out_plot_name=out_files["knee_plot"]
            )
            if min_norm_gini is None:
                min_norm_gini, min_norm_entropy = (
                    args.min_norm_gini,
                    args.min_norm_entropy,
                )
        else:
            min_norm_gini = 0.5
            min_norm_entropy = 0.75
            logging.warning(
                f"There's only one genome. Using min_norm_gini={min_norm_gini} and min_norm_entropy={min_norm_entropy}. Please check the results."
            )
    else:
        min_norm_gini, min_norm_entropy = args.min_norm_gini, args.min_norm_entropy

    filter_conditions = {
        "min_read_length": args.min_read_length,
        "min_read_count": args.min_read_count,
        "min_expected_breadth_ratio": args.min_expected_breadth_ratio,
        "min_breadth": args.min_breadth,
        "min_read_ani": args.min_read_ani,
        "min_coverage_evenness": args.min_coverage_evenness,
        "min_norm_entropy": min_norm_entropy,
        "min_norm_gini": min_norm_gini,
    }
    if args.only_stats:
        logging.info("Skipping filtering...")
    else:
        filter_reference_BAM(
            bam=args.bam,
            df=data_df,
            filter_conditions=filter_conditions,
            threads=args.threads,
            out_files=out_files,
            only_stats_filtered=args.only_stats_filtered,
            sort_memory=args.sort_memory,
            sort_by_name=args.sort_by_name,
        )
    logging.info(f"ALL DONE.")


if __name__ == "__main__":
    main()
