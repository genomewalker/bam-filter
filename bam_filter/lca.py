"""
CLI wrapper for the Cython-based LCA implementation.

The heavy lifting lives in :mod:`bam_filter.processor_lca`; this module simply
adapts parsed CLI arguments to the new engine while keeping backwards
compatible output handling.
"""

import gzip
from pathlib import Path
from collections import defaultdict

from bam_filter import logging as bf_logging
from bam_filter.processor_lca import run_lca
from bam_filter.lca_stats import run_lca_stats
from bam_filter.utils import create_output_files

LOG_TAG = "LCA"


def _info(message: str, *args) -> None:
    bf_logging.log(LOG_TAG, message, *args)


def _warn(message: str, *args) -> None:
    bf_logging.warn(message, *args)


def _classify_subtrees(
    summary_path: Path,
    taxdb,
    score_threshold: float = 0.0,
    pvalue_threshold: float = 0.05,
    authentic_fraction_threshold: float = 0.8,
    contaminant_fraction_threshold: float = 0.2,
) -> None:
    """
    Add subtree classification columns to LCA summary.

    For each taxon, computes both count-based and read-weighted metrics:
    - taxon_status: 'authentic' or 'contaminant' based on score/pvalue thresholds
    - subtree_status: 'authentic', 'contaminant', or 'mixed' based on read-weighted fraction
    - Count-based: authentic_descendant_count, authentic_descendant_fraction
    - Read-weighted: authentic_descendant_reads, total_descendant_reads, authentic_read_fraction

    Read weighting uses log10(1 + n_reads) to prevent single high-read taxa from dominating.

    Parameters
    ----------
    summary_path : Path
        LCA summary file (may be gzipped), modified in place
    taxdb : TaxonomyDatabase
        Loaded taxonomy database for parent lookups
    score_threshold : float
        Minimum authenticity_score for authentic classification (default 0.0)
    pvalue_threshold : float
        Maximum authenticity_pvalue for authentic classification (default 0.05)
    authentic_fraction_threshold : float
        Minimum read-weighted fraction for 'authentic' subtree status (default 0.8)
    contaminant_fraction_threshold : float
        Maximum read-weighted fraction for 'contaminant' subtree status (default 0.2)
    """
    import math

    is_gzipped = summary_path.name.endswith(".gz")
    open_func = gzip.open if is_gzipped else open

    # First pass: read all rows and classify each taxon
    rows = []
    taxid_to_idx = {}
    taxon_status = {}  # taxid -> 'authentic' or 'contaminant'
    taxon_reads = {}   # taxid -> n_reads
    taxon_weight = {}  # taxid -> log10(1 + n_reads)

    with open_func(summary_path, "rt") as fin:
        header = fin.readline()
        if not header:
            _warn("Empty LCA summary file: %s", summary_path)
            return

        columns = header.rstrip("\n").split("\t")

        try:
            taxid_idx = columns.index("taxid")
            score_idx = columns.index("authenticity_score")
            pvalue_idx = columns.index("authenticity_pvalue")
            reads_idx = columns.index("n_reads")
        except ValueError as e:
            _warn("Missing required column for classification: %s", e)
            return

        for line in fin:
            fields = line.rstrip("\n").split("\t")
            rows.append(fields)

            taxid = int(fields[taxid_idx])
            taxid_to_idx[taxid] = len(rows) - 1

            try:
                n_reads = int(fields[reads_idx])
            except (ValueError, IndexError):
                n_reads = 0
            taxon_reads[taxid] = n_reads
            taxon_weight[taxid] = math.log10(1 + n_reads)

            try:
                score = float(fields[score_idx])
                pvalue = float(fields[pvalue_idx])
            except (ValueError, IndexError):
                taxon_status[taxid] = "contaminant"
                continue

            if score >= score_threshold and pvalue <= pvalue_threshold:
                taxon_status[taxid] = "authentic"
            else:
                taxon_status[taxid] = "contaminant"

    # Build parent->children mapping from taxdb
    children = defaultdict(list)
    for taxid in taxid_to_idx:
        parent_id = taxdb.get_parent(taxid)
        if parent_id and parent_id != taxid and parent_id in taxid_to_idx:
            children[parent_id].append(taxid)

    # Compute subtree statistics for each taxon
    # Returns (auth_count, total_count, auth_reads, total_reads, auth_weight, total_weight)
    subtree_cache = {}

    def get_subtree_stats(taxid):
        if taxid in subtree_cache:
            return subtree_cache[taxid]

        is_auth = taxon_status.get(taxid) == "authentic"
        n_reads = taxon_reads.get(taxid, 0)
        weight = taxon_weight.get(taxid, 0.0)

        auth_count = 1 if is_auth else 0
        total_count = 1
        auth_reads = n_reads if is_auth else 0
        total_reads = n_reads
        auth_weight = weight if is_auth else 0.0
        total_weight = weight

        for child_id in children.get(taxid, []):
            c_auth_count, c_total_count, c_auth_reads, c_total_reads, c_auth_weight, c_total_weight = get_subtree_stats(child_id)
            auth_count += c_auth_count
            total_count += c_total_count
            auth_reads += c_auth_reads
            total_reads += c_total_reads
            auth_weight += c_auth_weight
            total_weight += c_total_weight

        subtree_cache[taxid] = (auth_count, total_count, auth_reads, total_reads, auth_weight, total_weight)
        return subtree_cache[taxid]

    for taxid in taxid_to_idx:
        get_subtree_stats(taxid)

    # Add new columns to header
    new_columns = [
        "taxon_status",
        "subtree_status",
        "authentic_descendant_count",
        "total_descendant_count",
        "authentic_descendant_fraction",
        "authentic_descendant_reads",
        "total_descendant_reads",
        "authentic_read_fraction",
    ]
    new_header = header.rstrip("\n") + "\t" + "\t".join(new_columns) + "\n"

    # Write updated file
    with open_func(summary_path, "wt") as fout:
        fout.write(new_header)

        for fields in rows:
            taxid = int(fields[taxid_idx])
            auth_count, total_count, auth_reads, total_reads, auth_weight, total_weight = subtree_cache.get(
                taxid, (0, 1, 0, 0, 0.0, 0.0)
            )

            # Descendant stats exclude self
            is_auth = taxon_status.get(taxid) == "authentic"
            self_reads = taxon_reads.get(taxid, 0)
            self_weight = taxon_weight.get(taxid, 0.0)

            desc_auth_count = auth_count - (1 if is_auth else 0)
            desc_total_count = total_count - 1
            desc_auth_reads = auth_reads - (self_reads if is_auth else 0)
            desc_total_reads = total_reads - self_reads
            desc_auth_weight = auth_weight - (self_weight if is_auth else 0.0)
            desc_total_weight = total_weight - self_weight

            # Count-based fraction
            if desc_total_count > 0:
                desc_count_fraction = desc_auth_count / desc_total_count
            else:
                desc_count_fraction = 0.0

            # Read-weighted fraction (using log weights)
            if desc_total_weight > 0:
                desc_read_fraction = desc_auth_weight / desc_total_weight
            else:
                desc_read_fraction = 0.0

            # Subtree status based on read-weighted fraction
            if desc_total_count == 0:
                # Leaf node: status equals taxon status
                subtree_stat = taxon_status.get(taxid, "contaminant")
            elif desc_read_fraction >= authentic_fraction_threshold:
                subtree_stat = "authentic"
            elif desc_read_fraction <= contaminant_fraction_threshold:
                subtree_stat = "contaminant"
            else:
                subtree_stat = "mixed"

            new_fields = [
                taxon_status.get(taxid, "contaminant"),
                subtree_stat,
                str(desc_auth_count),
                str(desc_total_count),
                f"{desc_count_fraction:.4f}",
                str(desc_auth_reads),
                str(desc_total_reads),
                f"{desc_read_fraction:.4f}",
            ]
            fout.write("\t".join(fields) + "\t" + "\t".join(new_fields) + "\n")

    # Summary statistics
    n_authentic = sum(1 for s in taxon_status.values() if s == "authentic")
    n_contaminant = sum(1 for s in taxon_status.values() if s == "contaminant")
    _info(
        "Subtree classification: %d authentic, %d contaminant taxa "
        "(score >= %.2f, pvalue <= %.2f)",
        n_authentic, n_contaminant, score_threshold, pvalue_threshold
    )


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

        # Add subtree classification columns
        score_threshold = getattr(args, "authenticity_score_threshold", 0.0)
        pvalue_threshold = getattr(args, "authenticity_pvalue_threshold", 0.05)
        subtree_auth_threshold = getattr(args, "subtree_authentic_threshold", 0.8)
        subtree_cont_threshold = getattr(args, "subtree_contaminant_threshold", 0.2)
        _classify_subtrees(
            summary_path,
            taxdb,
            score_threshold=score_threshold,
            pvalue_threshold=pvalue_threshold,
            authentic_fraction_threshold=subtree_auth_threshold,
            contaminant_fraction_threshold=subtree_cont_threshold,
        )

        if temp_per_read_path:
            try:
                Path(temp_per_read_path).unlink(missing_ok=True)
            except Exception:
                _warn(
                    "Failed to remove temporary per-read file %s",
                    temp_per_read_path,
                )
