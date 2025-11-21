"""High-level wrapper for BAM to Parquet conversion.

This module keeps a small Python façade around the optimized Cython
implementation in :mod:`bam_filter.processor_parquet_writer`. It preserves the
public API that other tools expect while delegating all heavy work to the
compiled extension (no pysam dependency).
"""

from pathlib import Path
from typing import Dict

import pyarrow as pa
import time

from bam_filter import logging as bf_logging
from bam_filter.processor_parquet_writer import (
    convert_bam_to_parquet as _cy_convert_bam_to_parquet,
    create_references_table as _cy_create_references_table,
)

LOG_TAG = "BAM-TO-PARQUET"


def _info(message: str, *args) -> None:
    bf_logging.log(LOG_TAG, message, *args)


def _warn(message: str, *args) -> None:
    bf_logging.warn(message, *args)


def _announce_stage(title: str, detail: str = "") -> None:
    bf_logging.summary("")
    bf_logging.summary("┌─ %s", title)
    if detail:
        bf_logging.summary("│ %s", detail)
    bf_logging.summary("└─────────────────────────────────────────────────────────────")


def _log_stage(stage: str, start_time: float, level: int = 0) -> None:
    duration = time.perf_counter() - start_time
    bf_logging.verbose(level, LOG_TAG, "stage=%s duration=%.2fs", stage, duration)


def get_parquet_schema(
    include_read_names: bool = False, include_sequences: bool = True
) -> pa.Schema:
    """Return the PyArrow schema used for alignment Parquet files."""
    fields = [
        ("read_id", pa.uint32()),
        ("ref_id", pa.uint32()),
        ("position", pa.int32()),
        ("end_position", pa.int32()),
        ("mapq", pa.uint8()),
        ("flag", pa.uint16()),
        ("ani", pa.float32()),
        ("alignment_score", pa.float32()),
        ("pmd_score", pa.float32()),
        ("num_mismatches", pa.uint32()),
        ("alignment_length", pa.uint32()),
        ("template_length", pa.int32()),
        ("mate_ref_id", pa.int32()),
        ("mate_position", pa.int32()),
    ]

    if include_read_names:
        fields.append(("read_name", pa.string()))

    fields.append(("cigar", pa.string()))

    if include_sequences:
        fields.append(("sequence", pa.string()))

    fields.append(("quality", pa.binary()))
    fields.append(("tags", pa.binary()))

    return pa.schema(fields)


def get_references_schema() -> pa.Schema:
    """Return the PyArrow schema for the reference dimension table."""
    return pa.schema(
        [
            ("ref_id", pa.uint32()),
            ("ref_name", pa.string()),
            ("ref_length", pa.uint32()),
            ("ref_partition", pa.uint16()),
        ]
    )


def create_references_table(bam_path: str, num_partitions: int = 256) -> pa.Table:
    """Build a references table from the BAM header using the Cython reader."""
    _info("Creating references table from %s", bam_path)
    table = _cy_create_references_table(bam_path, num_partitions)
    _info("Created references table with %d references", table.num_rows)
    return table


def convert_bam_to_parquet(
    bam_path: str,
    output_base_path: str,
    num_partitions: int = -1,
    batch_size: int = -1,
    num_threads: int = 1,
    compression: str = "zstd",
    compression_level: int = 3,
    include_read_names: bool = True,
    include_sequences: bool = True,
    include_sequence_text: bool = False,
    calculate_pmd: bool = True,
    min_read_length: int = 0,
    max_read_length: int = 0,
    min_read_ani: float = 0.0,
    min_mapq: int = 0,
) -> Dict[str, int]:
    """Convert BAM to partitioned Parquet format (no filtering applied).

    Parameters other than ``calculate_pmd`` are preserved for backward
    compatibility but do not influence filtering—the Cython converter emits all
    mapped alignments.
    """
    if not include_sequences:
        _warn("Sequences are always included in Parquet output; overriding --no-sequences")
        include_sequences = True

    output_path = Path(output_base_path)
    _info("Starting BAM to Parquet conversion: %s", bam_path)
    _info("Output directory: %s", output_path.resolve())
    auto_partitions = num_partitions <= 0
    auto_batch = batch_size <= 0
    partition_desc = "auto" if auto_partitions else str(num_partitions)
    batch_desc = "auto" if auto_batch else str(batch_size)
    _info("Partitions: %s, Batch size: %s", partition_desc, batch_desc)
    _info(
        "Compression: %s (level %d) | Read names: %s | Sequences: %s | Sequence text: %s | PMD: %s",
        compression,
        compression_level,
        include_read_names,
        include_sequences,
        include_sequence_text,
        calculate_pmd,
    )

    if (
        min_read_length
        or max_read_length
        or min_read_ani
        or min_mapq
    ):
        _warn("Read-level filters are ignored during Parquet conversion; all alignments are emitted.")

    _announce_stage("Setup", "Preparing output directories and schema files")
    stage_timer = time.perf_counter()
    output_path.mkdir(parents=True, exist_ok=True)

    alignments_dir = output_path / "alignments"
    alignments_dir.mkdir(exist_ok=True)

    references_dir = output_path / "references"
    references_dir.mkdir(exist_ok=True)
    _log_stage("setup", stage_timer)

    _announce_stage("Streaming", "Reading BAM alignments and writing Parquet partitions")
    stage_timer = time.perf_counter()
    stats = _cy_convert_bam_to_parquet(
        bam_path=bam_path,
        output_base_path=str(output_path),
        num_partitions=num_partitions,
        batch_size=batch_size,
        num_threads=num_threads,
        compression=compression,
        compression_level=compression_level,
        include_read_names=include_read_names,
        include_sequences=include_sequences,
        include_sequence_text=include_sequence_text,
        calculate_pmd=calculate_pmd,
    )
    _log_stage("streaming", stage_timer)

    _announce_stage("Summary", "Recording conversion statistics and metadata")
    stage_timer = time.perf_counter()
    if stats.get("auto_num_partitions", False):
        _info("Auto-selected partitions: %s", stats.get("num_partitions_used"))
    else:
        _info("Partitions used: %s", stats.get("num_partitions_used"))
    if stats.get("auto_batch_size", False):
        _info("Auto-selected batch size: %s", stats.get("batch_size_used"))
    else:
        _info("Batch size used: %s", stats.get("batch_size_used"))
    if "estimated_total_alignments" in stats:
        _info(
            "Estimated total alignments from index: %s",
            stats["estimated_total_alignments"],
        )
    _info("Conversion complete!")
    _info("  Total alignments processed: %d", stats.get("total_alignments", 0))
    _info("  Alignments written: %d", stats.get("written_alignments", 0))
    _info("  Partitions created: %d", stats.get("partitions_created", 0))
    _log_stage("summary", stage_timer)

    return stats


def do_bam_to_parquet(args) -> Dict[str, int]:
    """CLI entry point wrapper used by ``filterBAM``."""
    return convert_bam_to_parquet(
        bam_path=args.bam,
        output_base_path=args.output,
        num_partitions=getattr(args, "num_partitions", -1),
        batch_size=getattr(args, "batch_size", -1),
        num_threads=getattr(args, "threads", 1),
        compression=getattr(args, "compression", "zstd"),
        compression_level=getattr(args, "compression_level", 3),
        include_read_names=getattr(args, "include_read_names", True),
        include_sequences=getattr(args, "include_sequences", True),
        include_sequence_text=getattr(args, "include_sequence_text", False),
        calculate_pmd=getattr(args, "calculate_pmd", True),
        min_read_length=getattr(args, "min_read_length", 0),
        max_read_length=getattr(args, "max_read_length", 0),
        min_read_ani=getattr(args, "min_read_ani", 0.0),
        min_mapq=getattr(args, "min_mapq", 0),
    )
