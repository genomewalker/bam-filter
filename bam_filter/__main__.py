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

from bam_filter.sam_utils import processBam, filterReferenceBAMfile

from bam_filter.utils import (
    get_arguments,
)

log = logging.getLogger("my_logger")


def main():

    logging.basicConfig(
        level=logging.DEBUG, format="%(levelname)s ::: %(asctime)s ::: %(message)s"
    )

    args = get_arguments()
    logging.getLogger("my_logger").setLevel(
        logging.DEBUG if args.debug else logging.INFO
    )
    datadf = processBam(bam=args.bam, threads=args.threads)
    datadf_filtered = datadf.loc[
        (datadf["breadthExpRatio"] > args.breadthExpRatio)
        & (datadf["covEvenness"] > args.covEvenness)
        & (datadf["nReads"] > args.nReads)
    ]
    datadf.to_csv(args.output, sep="\t", index=False)
    refs_dict = dict(zip(datadf["chromosome"], datadf["referenceLength"]))
    filterReferenceBAMfile(args.bam, refs_dict, outBAMfile="test.bam")


if __name__ == "__main__":
    main()
