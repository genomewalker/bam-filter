#!/usr/bin/env python3
"""
Fast Parquet converter using PyArrow directly.

Key optimizations:
- Uses PyArrow's efficient Table API
- Batch processing (100K records at a time)
- Direct Parquet writing (no DuckDB overhead)
- Minimal Python overhead with numpy arrays
"""

import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path
import time


# Schema matching DuckDB version
SCHEMA = pa.schema([
    ('read_id', pa.uint64()),
    ('read_name', pa.string()),
    ('ref_id', pa.uint32()),
    ('ref_name', pa.string()),
    ('position', pa.int32()),
    ('end_position', pa.int32()),
    ('mapq', pa.uint8()),
    ('flag', pa.uint16()),
    ('alignment_length', pa.uint16()),
    ('template_length', pa.int32()),
    ('mate_ref_id', pa.int32()),
    ('mate_position', pa.int32()),
    ('alignment_score', pa.int32()),
    ('xs_score', pa.int32()),
    ('edit_distance', pa.uint16()),
    ('num_mismatches', pa.uint8()),
    ('num_gap_opens', pa.uint8()),
    ('num_gap_extensions', pa.uint8()),
    ('md_string', pa.string()),
    ('ani', pa.float32()),
    ('pmd_score', pa.float32()),
    ('zs_score', pa.float32()),
    ('zp_posterior', pa.float32()),
    ('lca_taxid', pa.int32()),
    ('reassigned_ref_id', pa.int32()),
    ('filter_passed', pa.bool_()),
    ('read_group', pa.string()),
    ('cigar', pa.string()),
    ('sequence', pa.binary()),
    ('quality', pa.binary()),
])


def convert_sam_bam_to_parquet_arrow(
    input_file: str,
    output_dir: str,
    num_partitions: int = 128,
    batch_size: int = 100000,
    write_by_reference: bool = True,
    write_by_read: bool = True,
):
    """
    Convert SAM/BAM to Parquet using PyArrow directly.

    Much faster than DuckDB approach because:
    - No database file creation overhead
    - Direct columnar writes
    - Optimized Arrow C++ backend
    """
    import pysam
    from collections import defaultdict

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    if write_by_reference:
        (output_path / "by_reference").mkdir(exist_ok=True)
    if write_by_read:
        (output_path / "by_read").mkdir(exist_ok=True)

    # Partition writers
    ref_writers = {}
    read_writers = {}

    # Partition buffers (lists of record dicts)
    ref_buffers = defaultdict(list)
    read_buffers = defaultdict(list)

    start_time = time.time()
    total_records = 0

    # Open SAM/BAM file
    samfile = pysam.AlignmentFile(input_file, "r")

    print(f"Processing {input_file}...")
    last_print = time.time()

    for aln in samfile:
        if aln.reference_id < 0:
            continue

        # Extract record data
        record = {
            'read_id': total_records,
            'read_name': aln.query_name,
            'ref_id': aln.reference_id,
            'ref_name': aln.reference_name,
            'position': aln.reference_start,
            'end_position': aln.reference_end,
            'mapq': aln.mapping_quality,
            'flag': aln.flag,
            'alignment_length': aln.query_length or 0,
            'template_length': aln.template_length,
            'mate_ref_id': aln.next_reference_id if aln.next_reference_id >= 0 else -1,
            'mate_position': aln.next_reference_start if aln.next_reference_start >= 0 else -1,
            'alignment_score': aln.get_tag('AS') if aln.has_tag('AS') else None,
            'xs_score': aln.get_tag('XS') if aln.has_tag('XS') else None,
            'edit_distance': aln.get_tag('NM') if aln.has_tag('NM') else None,
            'num_mismatches': aln.get_tag('XM') if aln.has_tag('XM') else None,
            'num_gap_opens': aln.get_tag('XO') if aln.has_tag('XO') else None,
            'num_gap_extensions': aln.get_tag('XG') if aln.has_tag('XG') else None,
            'md_string': aln.get_tag('MD') if aln.has_tag('MD') else None,
            'ani': None,
            'pmd_score': aln.get_tag('PMD') if aln.has_tag('PMD') else (aln.get_tag('PM') if aln.has_tag('PM') else None),
            'zs_score': aln.get_tag('ZS') if aln.has_tag('ZS') else None,
            'zp_posterior': aln.get_tag('ZP') if aln.has_tag('ZP') else None,
            'lca_taxid': aln.get_tag('ZT') if aln.has_tag('ZT') else None,
            'reassigned_ref_id': aln.get_tag('ZR') if aln.has_tag('ZR') else None,
            'filter_passed': True,
            'read_group': aln.get_tag('RG') if aln.has_tag('RG') else None,
            'cigar': aln.cigarstring,
            'sequence': aln.query_sequence.encode() if aln.query_sequence else b'',
            'quality': bytes(aln.query_qualities) if aln.query_qualities is not None else b'',
        }

        # Calculate ANI
        if record['edit_distance'] is not None and record['alignment_length'] > 0:
            record['ani'] = (1.0 - record['edit_distance'] / record['alignment_length']) * 100.0

        # Add to partition buffers
        ref_partition = aln.reference_id % num_partitions
        read_partition = total_records % num_partitions

        if write_by_reference:
            ref_buffers[ref_partition].append(record)

            # Flush if buffer full
            if len(ref_buffers[ref_partition]) >= batch_size:
                _flush_partition(ref_buffers[ref_partition], ref_partition, ref_writers,
                               output_path / "by_reference", True)
                ref_buffers[ref_partition] = []

        if write_by_read:
            read_buffers[read_partition].append(record)

            if len(read_buffers[read_partition]) >= batch_size:
                _flush_partition(read_buffers[read_partition], read_partition, read_writers,
                               output_path / "by_read", False)
                read_buffers[read_partition] = []

        total_records += 1

        # Progress
        if time.time() - last_print >= 10:
            elapsed = time.time() - start_time
            print(f"\rProgress: {total_records:,} records ({total_records/elapsed/1e6:.2f} M/sec, {elapsed:.0f}s)",
                  end='', flush=True)
            last_print = time.time()

    samfile.close()

    # Flush remaining buffers
    print(f"\n\nFlushing remaining buffers...")
    for partition, buffer in ref_buffers.items():
        if buffer:
            _flush_partition(buffer, partition, ref_writers, output_path / "by_reference", True)

    for partition, buffer in read_buffers.items():
        if buffer:
            _flush_partition(buffer, partition, read_writers, output_path / "by_read", False)

    # Close all writers
    for writer in ref_writers.values():
        writer.close()
    for writer in read_writers.values():
        writer.close()

    duration = time.time() - start_time

    return {
        'total_records': total_records,
        'processing_time_seconds': duration,
    }


def _flush_partition(buffer, partition_id, writers, output_dir, is_reference):
    """Flush buffer to Parquet file."""
    if not buffer:
        return

    # Convert list of dicts to columnar format
    columns = {field.name: [] for field in SCHEMA}

    for record in buffer:
        for field in SCHEMA:
            columns[field.name].append(record[field.name])

    # Create Arrow arrays
    arrays = []
    for field in SCHEMA:
        col_data = columns[field.name]
        arrays.append(pa.array(col_data, type=field.type))

    # Create table
    table = pa.Table.from_arrays(arrays, schema=SCHEMA)

    # Get or create writer
    if partition_id not in writers:
        prefix = "by_reference" if is_reference else "by_read"
        filename = output_dir / f"{prefix}_p{partition_id:04d}.parquet"
        writers[partition_id] = pq.ParquetWriter(
            filename,
            SCHEMA,
            compression='zstd',
            use_dictionary=True,
        )

    # Write table
    writers[partition_id].write_table(table)
