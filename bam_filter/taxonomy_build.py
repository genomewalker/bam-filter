"""
Taxonomy build utilities and CLI entry points.

This module provides the `build_taxonomy` command used by the CLI. It mirrors
the naming convention adopted by other subsystems (for example, `processor_*`
and `stats_*`) so taxonomy-related helpers live under a consistent prefix.

Users can build a Parquet-backed taxonomy database directly from the NCBI
taxdump files, and convert accession→taxid tables on the fly without relying
on additional standalone scripts.
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Iterable, Optional, Tuple

from bam_filter import logging as bf_logging
from bam_filter.taxonomy import TaxonomyDB


def _lazy_import_pyarrow():
    """
    Import pyarrow lazily so the CLI only pays the import cost when needed.
    """

    try:
        import pyarrow as pa  # type: ignore
        import pyarrow.csv as csv  # type: ignore
        import pyarrow.parquet as pq  # type: ignore
    except ModuleNotFoundError as exc:  # pragma: no cover - defensive import guard
        raise RuntimeError(
            "pyarrow is required for taxonomy Parquet conversion. "
            "Install it via `pip install pyarrow` or `conda install pyarrow`."
        ) from exc
    return pa, csv, pq


def convert_gz_to_parquet(
    gz_file: str,
    parquet_file: str,
    batch_size: int = 10_000_000,
    *,
    compression: Optional[str] = "zstd",
    compression_level: Optional[int] = 3,
) -> int:
    """
    Convert an acc2taxid `.gz` file to Parquet using Arrow's streaming CSV reader.

    Parameters
    ----------
    gz_file : str
        Path to the gzipped accession→taxid TSV file.
    parquet_file : str
        Output Parquet file.
    batch_size : int, optional
        Number of rows per Parquet row group (default: 10M).
    compression : str, optional
        Compression codec to use (default: zstd).
    compression_level : int, optional
        Compression level passed to the Parquet writer.

    Returns
    -------
    int
        Total number of rows written.
    """

    gz_path = Path(gz_file)
    parquet_path = Path(parquet_file)

    bf_logging.log("BUILD-TAXONOMY", "Converting %s → %s", gz_path.name, parquet_path.name)
    start_time = time.time()

    pa, csv, pq = _lazy_import_pyarrow()

    with gzip.open(gz_path, "rt") as handle:
        header = handle.readline().strip()
        columns = header.split("\t")

    if "accession.version" in columns:
        acc_column = "accession.version"
    elif "accession" in columns:
        acc_column = "accession"
    else:
        raise RuntimeError(
            f"Could not find an accession column in {gz_path}. "
            "Expected 'accession.version' or 'accession'."
        )

    if "taxid" not in columns:
        raise RuntimeError(f"Could not find a 'taxid' column in {gz_path}.")

    decompress_start = time.time()
    with tempfile.NamedTemporaryFile(mode="wb", suffix=".tsv", delete=False) as tmp_file:
        tmp_path = Path(tmp_file.name)

    pigz_path = shutil.which("pigz")
    try:
        if pigz_path:
            with open(tmp_path, "wb") as dest:
                subprocess.run([pigz_path, "-dc", str(gz_path)], stdout=dest, check=True)
        else:
            with gzip.open(gz_path, "rb") as source, open(tmp_path, "wb") as dest:
                shutil.copyfileobj(source, dest, length=10 * 1024 * 1024)
    except (OSError, subprocess.CalledProcessError) as exc:  # pragma: no cover - defensive
        tmp_path.unlink(missing_ok=True)  # type: ignore[attr-defined]
        raise RuntimeError(f"Failed to decompress {gz_path}: {exc}") from exc

    bf_logging.log(
        "BUILD-TAXONOMY",
        "  Decompressed in %.1fs using %s",
        time.time() - decompress_start,
        "pigz" if pigz_path else "Python gzip",
    )

    read_options = csv.ReadOptions(use_threads=True, block_size=32 * 1024 * 1024)
    parse_options = csv.ParseOptions(delimiter="\t")
    convert_options = csv.ConvertOptions(
        include_columns=[acc_column, "taxid"],
        column_types={acc_column: pa.string(), "taxid": pa.int32()},
    )

    parquet_path.parent.mkdir(parents=True, exist_ok=True)
    writer = None
    total_rows = 0
    try:
        reader = csv.open_csv(
            tmp_path,
            read_options=read_options,
            parse_options=parse_options,
            convert_options=convert_options,
        )

        for batch in reader:
            table = pa.Table.from_batches([batch]).rename_columns(["accession", "taxid"])
            if writer is None:
                writer = pq.ParquetWriter(
                    parquet_path,
                    table.schema,
                    compression=compression,
                    compression_level=compression_level,
                    use_dictionary=True,
                )
            writer.write_table(table, row_group_size=batch_size)
            total_rows += batch.num_rows
    finally:
        if writer is not None:
            writer.close()
        tmp_path.unlink(missing_ok=True)  # type: ignore[attr-defined]

    elapsed = time.time() - start_time
    bf_logging.log(
        "BUILD-TAXONOMY",
        f"  ↳ {total_rows:,d} rows written in {elapsed:.1f}s",
    )
    return total_rows


def _merge_parquet_tables(parquet_files: Iterable[str], output_path: str) -> int:
    """
    Merge multiple Parquet accession maps into a single file.
    """

    pa, _, pq = _lazy_import_pyarrow()
    tables = []
    for parquet_file in parquet_files:
        bf_logging.log("BUILD-TAXONOMY", "  Reading %s...", os.path.basename(parquet_file))
        tables.append(pq.read_table(parquet_file))

    merged = pa.concat_tables(tables, promote=True)
    pq.write_table(
        merged,
        output_path,
        compression="zstd",
        compression_level=3,
        use_dictionary=True,
    )
    return merged.num_rows


def _count_parquet_rows(parquet_file: str) -> int:
    """
    Count rows in a Parquet file without bringing the entire table into memory.
    """

    if not os.path.exists(parquet_file):
        return 0
    _, _, pq = _lazy_import_pyarrow()
    return pq.read_table(parquet_file, columns=["accession"]).num_rows


def build_taxonomy(args) -> None:
    """
    Build taxonomy database from NCBI dump files and save accompanying Parquet data.
    """

    bf_logging.log("BUILD-TAXONOMY", "Starting taxonomy database construction")
    bf_logging.log("BUILD-TAXONOMY", "")

    if not os.path.exists(args.nodes):
        bf_logging.error("BUILD-TAXONOMY: nodes file not found: %s", args.nodes)
        sys.exit(1)
    if not os.path.exists(args.names):
        bf_logging.error("BUILD-TAXONOMY: names file not found: %s", args.names)
        sys.exit(1)

    for acc2taxid_file in args.acc2taxid:
        if not os.path.exists(acc2taxid_file):
            bf_logging.error("BUILD-TAXONOMY: acc2taxid file not found: %s", acc2taxid_file)
            sys.exit(1)

    os.makedirs(args.output, exist_ok=True)

    bf_logging.log("BUILD-TAXONOMY", "Step 1/3: Loading taxonomy tree (nodes + names)")
    bf_logging.log("BUILD-TAXONOMY", "  nodes: %s", args.nodes)
    bf_logging.log("BUILD-TAXONOMY", "  names: %s", args.names)

    try:
        tax = TaxonomyDB.from_ncbi(
            nodes_file=args.nodes,
            names_file=args.names,
            acc2taxid_file=None,
            num_threads=args.num_threads,
        )
        bf_logging.log("BUILD-TAXONOMY", f"  ✓ Loaded {tax.n_nodes:,d} taxonomy nodes")
    except Exception as exc:  # pragma: no cover - defensive
        bf_logging.error("BUILD-TAXONOMY: Error loading taxonomy: %s", exc)
        sys.exit(1)

    bf_logging.log("BUILD-TAXONOMY", "")
    bf_logging.log(
        "BUILD-TAXONOMY",
        "Step 2/3: Processing accession mappings (%d file(s))",
        len(args.acc2taxid),
    )

    parquet_files = []
    for index, acc2taxid_file in enumerate(args.acc2taxid):
        bf_logging.log(
            "BUILD-TAXONOMY",
            "  File %d/%d: %s",
            index + 1,
            len(args.acc2taxid),
            acc2taxid_file,
        )

        if len(args.acc2taxid) == 1:
            parquet_file = os.path.join(args.output, "accession_map.parquet")
        else:
            base_name = os.path.basename(acc2taxid_file).replace(".gz", ".parquet")
            parquet_file = os.path.join(args.output, f"accession_map_{index:02d}_{base_name}")

        if acc2taxid_file.endswith(".parquet"):
            bf_logging.log("BUILD-TAXONOMY", "    Copying existing Parquet file…")
            shutil.copy2(acc2taxid_file, parquet_file)
            parquet_files.append(parquet_file)
        elif acc2taxid_file.endswith(".gz"):
            if os.path.exists(parquet_file):
                bf_logging.log("BUILD-TAXONOMY", "    Using cached Parquet file: %s", parquet_file)
            else:
                try:
                    convert_gz_to_parquet(acc2taxid_file, parquet_file)
                except Exception as exc:  # pragma: no cover - defensive
                    bf_logging.error("BUILD-TAXONOMY: Error converting %s: %s", acc2taxid_file, exc)
                    sys.exit(1)
            parquet_files.append(parquet_file)
        else:
            bf_logging.error("BUILD-TAXONOMY: Unsupported file format: %s", acc2taxid_file)
            bf_logging.error("BUILD-TAXONOMY: Only .parquet and .gz files are supported")
            sys.exit(1)

    bf_logging.log("BUILD-TAXONOMY", "")

    if args.cache_taxids:
        bf_logging.log(
            "BUILD-TAXONOMY",
            f"Building LCA cache for {len(args.cache_taxids):,d} taxids",
        )
        try:
            tax.build_lca_cache(args.cache_taxids)
        except Exception as exc:  # pragma: no cover - defensive
            bf_logging.warn("BUILD-TAXONOMY: Failed to build LCA cache: %s", exc)

    bf_logging.log("BUILD-TAXONOMY", "Step 3/3: Saving taxonomy database")
    try:
        tax.save(args.output)

        if len(parquet_files) > 1:
            bf_logging.log(
                "BUILD-TAXONOMY",
                "Merging %d accession Parquet shards…",
                len(parquet_files),
            )
            final_parquet = os.path.join(args.output, "accession_map.parquet")
            merged_rows = _merge_parquet_tables(
                parquet_files,
                final_parquet,
            )
            bf_logging.log(
                "BUILD-TAXONOMY",
                f"  ✓ Merged {merged_rows:,d} accession mappings",
            )
            for temp_file in parquet_files:
                if os.path.abspath(temp_file) != os.path.abspath(final_parquet):
                    try:
                        os.remove(temp_file)
                    except OSError:
                        pass
    except Exception as exc:  # pragma: no cover - defensive
        bf_logging.error("BUILD-TAXONOMY: Error saving taxonomy: %s", exc)
        sys.exit(1)

    bf_logging.log("BUILD-TAXONOMY", "")
    bf_logging.log("BUILD-TAXONOMY", "✓ Taxonomy database built successfully!")
    bf_logging.log("BUILD-TAXONOMY", f"  Taxonomy nodes: {tax.n_nodes:,d}")

    accession_rows = _count_parquet_rows(os.path.join(args.output, "accession_map.parquet"))
    if accession_rows:
        bf_logging.log("BUILD-TAXONOMY", f"  Accession mappings: {accession_rows:,d}")

    if getattr(tax, "lca_cache", None) is not None:
        bf_logging.log(
            "BUILD-TAXONOMY",
            f"  LCA cache: {getattr(tax.lca_cache, 'n_cached', 0):,d} taxids",
        )

    bf_logging.log("BUILD-TAXONOMY", "")
    bf_logging.log("BUILD-TAXONOMY", "Database saved to: %s/", args.output)
    bf_logging.log("BUILD-TAXONOMY", "")
    bf_logging.log("BUILD-TAXONOMY", "To use this database:")
    bf_logging.log("BUILD-TAXONOMY", "  filterBAM lca --taxonomy %s ...", args.output)
