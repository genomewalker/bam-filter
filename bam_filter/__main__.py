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
import numpy as np
from bam_filter.utils import get_arguments, create_output_files

log = logging.getLogger("my_logger")


def main():

    logging.basicConfig(
        level=logging.DEBUG, format="%(levelname)s ::: %(asctime)s ::: %(message)s"
    )

    args = get_arguments()
    logging.getLogger("my_logger").setLevel(
        logging.DEBUG if args.debug else logging.INFO
    )

    out_files = create_output_files(prefix=args.prefix, bam=args.bam)

    data = process_bam(
        bam=args.bam, threads=args.threads, reference_lengths=args.reference_lengths
    )

    data_df = pd.DataFrame([x.to_summary() for x in data])

    logging.info(f"Writing reference statistics to {out_files['stats']}")
    data_df.to_csv(out_files["stats"], sep="\t", index=False, compression="gzip")

    filter_conditions = {
        "min_read_length": args.min_read_length,
        "min_read_count": args.min_read_count,
        "min_expected_breadth_ratio": args.min_expected_breadth_ratio,
        "min_read_ani": args.min_read_ani,
        "min_coverage_evenness": args.min_coverage_evenness,
    }

    filter_reference_BAM(
        bam=args.bam,
        df=data_df,
        filter_conditions=filter_conditions,
        threads=args.threads,
        out_files=out_files,
        sort_memory=args.sort_memory,
    )
    logging.info(f"ALL DONE.")


if __name__ == "__main__":
    main()
