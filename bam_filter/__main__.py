#
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

from bam_filter.reassign import reassign
from bam_filter.filter import filter_references
from bam_filter.lca import do_lca
from bam_filter.build_taxonomy import build_taxonomy
from bam_filter.utils import (
    get_arguments,
)
from bam_filter import __version__
from bam_filter import logging as bf_logging


def main():
    args = get_arguments()
    bf_logging.set_level(getattr(args, "verbosity", bf_logging.get_level()))
    bf_logging.log("CLI", f"Using filterBAM version: {__version__}")
    if args.action == "reassign":
        reassign(args)
    elif args.action == "filter":
        filter_references(args)
    elif args.action == "lca":
        do_lca(args)
    elif args.action == "build-taxonomy":
        build_taxonomy(args)


if __name__ == "__main__":
    main()
