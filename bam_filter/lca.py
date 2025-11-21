"""
CLI wrapper for the Cython-based LCA implementation.

The heavy lifting lives in :mod:`bam_filter.processor_lca`; this module simply
adapts parsed CLI arguments to the new engine while keeping backwards
compatible output handling.
"""

from pathlib import Path

from bam_filter import logging as bf_logging
from bam_filter.processor_lca import run_lca
from bam_filter.lca_stats import run_lca_stats
from bam_filter.utils import create_output_files

LOG_TAG = "LCA"


def _info(message: str, *args) -> None:
    bf_logging.log(LOG_TAG, message, *args)


def _warn(message: str, *args) -> None:
    bf_logging.warn(message, *args)


def do_lca(args):
    """
    Execute the LCA subcommand using the optimised Cython pipeline.
    """
    prefix_arg = getattr(args, "prefix", None)
    if prefix_arg:
        prefix = prefix_arg
    else:
        prefix = Path(args.bam).with_suffix("").name

    out_files = create_output_files(
        bam=args.bam,
        prefix=prefix,
        tmp_dir=getattr(args, "tmp_dir", None),
        mode="lca",
        lca_summary=getattr(args, "lca_summary", None),
    )
    summary_path = Path(out_files["lca_summary"])
    tmp_dir_path = Path(out_files["tmp_dir"])

    taxonomy_db_path = Path(args.taxonomy_db)
    if not taxonomy_db_path.exists() or not taxonomy_db_path.is_dir():
        raise FileNotFoundError(
            f"Taxonomy database directory not found: {taxonomy_db_path}"
        )

    accession_map_path = taxonomy_db_path / "accession_map.parquet"
    if not accession_map_path.exists():
        raise FileNotFoundError(
            f"accession_map.parquet not found inside {taxonomy_db_path}. "
            "Run `filterBAM build-taxonomy` first or provide the accession map."
        )

    custom_flag = bool(getattr(args, "custom", False))
    if custom_flag:
        _warn(
            "--custom is deprecated; accession maps are expected in Parquet format. "
            "Continuing for backwards compatibility."
        )

    scale_arg = getattr(args, "scale", None)
    scale_factor = int(scale_arg) if scale_arg is not None else 1_000_000
    if scale_factor <= 0:
        raise ValueError("--scale must be a positive integer")

    _info("Input BAM: %s", args.bam)
    _info("Taxonomy database: %s", taxonomy_db_path)
    _info("Scale factor: %d", scale_factor)

    per_read_path = getattr(args, "lca_per_read", None)
    if per_read_path:
        _info("Per-read LCA output: %s", per_read_path)

    enable_stats = bool(getattr(args, "lca_stats", False))
    stats_tmp_path = None
    if enable_stats:
        if summary_path.name.endswith(".gz"):
            stats_tmp_name = summary_path.name[:-3] + ".tmp.gz"
        else:
            stats_tmp_name = f"{summary_path.name}.tmp"
        stats_tmp_path = summary_path.parent / stats_tmp_name

    temp_per_read_path = None
    if enable_stats and not per_read_path:
        temp_per_read_path = tmp_dir_path / f"{prefix}_per-read.tsv.gz"
        per_read_path = str(temp_per_read_path)
        _info(
            "Per-read LCA output not provided; writing to %s to support --stats",
            per_read_path,
        )

    if enable_stats:
        _info("LCA taxonomic stats will be merged into %s", summary_path)

    lca_result, taxdb = run_lca(
        bam_path=args.bam,
        output_path=str(summary_path),
        rank=getattr(args, "rank_lca", "genus"),
        custom_acc=custom_flag,
        threads=getattr(args, "threads", 1),
        min_read_ani=getattr(args, "min_read_ani", 0.0),
        min_read_length=getattr(args, "min_read_length", 30),
        min_read_count=getattr(args, "min_read_count", 1),
        scale=scale_factor,
        verbose=not getattr(args, "quiet", False),
        stats_path=None,
        reference_lengths_tsv=getattr(args, "reference_lengths_tsv", None),
        taxonomy_db_dir=str(taxonomy_db_path),
        per_read_path=per_read_path,
    )

    _info("LCA taxonomy summary written to %s", summary_path.resolve())
    if per_read_path:
        _info("Per-read LCA written to %s", Path(per_read_path).resolve())

    if enable_stats:
        print("")
        print("┌─ LCA Stats Phase: Computing per-taxon quality metrics")
        print("│ Calculating coverage, breadth, and quality statistics for each taxon")
        print("└─────────────────────────────────────────────────────────────")
        run_lca_stats(
            bam_path=args.bam,
            lca_per_read_path=per_read_path,
            output_path=str(stats_tmp_path),
            taxonomy_db_path=str(taxonomy_db_path),
            num_threads=getattr(args, "threads", 1),
            verbose=bf_logging.should_log(bf_logging.LogLevel.INFO),
            min_read_ani=getattr(args, "min_read_ani", 0.0),
            min_read_length=getattr(args, "min_read_length", 0),
            max_read_length=getattr(args, "max_read_length", (2**31) - 1),
            scale=scale_factor,
            trim_ends=int(getattr(args, "trim_ends", 0)),
            trim_min=getattr(args, "trim_min", 10),
            trim_max=getattr(args, "trim_max", 90),
            taxdb=taxdb,
        )
        Path(stats_tmp_path).replace(summary_path)
        bf_logging.summary("LCA stats merged into %s (added quality metrics)", summary_path.resolve())

        if temp_per_read_path:
            try:
                Path(temp_per_read_path).unlink(missing_ok=True)
            except Exception:
                _warn(
                    "Failed to remove temporary per-read file %s",
                    temp_per_read_path,
                )
