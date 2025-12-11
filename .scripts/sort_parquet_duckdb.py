#!/usr/bin/env python3
"""
Sort Parquet file by (ref_id, position) using DuckDB for optimal compression.

DuckDB can efficiently sort large files with minimal memory by using disk spilling.
Sorting clusters similar alignments together, improving:
1. RLE compression for ref_id (long runs per chromosome)
2. Delta compression for position (smaller deltas)
3. Dictionary encoding efficiency
"""

import os
import sys
import time
from pathlib import Path
from typing import Dict, Any, Optional

import duckdb


def sort_parquet_by_position(
    input_file: str,
    output_file: str,
    compression_level: int = 6,
    row_group_size: int = 500000,
    num_threads: int = 8,
    temp_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Sort a Parquet file by (ref_id, position) for better compression.

    Uses DuckDB's external merge sort for memory efficiency.

    Args:
        input_file: Input Parquet file
        output_file: Output sorted Parquet file
        compression_level: ZSTD compression level (1-22, default 6)
        row_group_size: Parquet row group size (default 500K)
        num_threads: Number of DuckDB threads
        temp_dir: Temporary directory for disk spilling

    Returns:
        Dictionary with timing and size statistics
    """
    print("=" * 70)
    print("DuckDB Parquet Position Sorter")
    print("=" * 70)
    print(f"Input:  {input_file}")
    print(f"Output: {output_file}")
    print(f"Compression level: {compression_level}")
    print(f"Row group size: {row_group_size:,}")
    print(f"Threads: {num_threads}")
    print()

    # Validate input
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")

    input_size = os.path.getsize(input_file)
    print(f"Input size: {input_size / 1024 / 1024:.1f} MB")

    # Create DuckDB connection with settings for large sort
    con = duckdb.connect(":memory:")
    con.execute(f"SET threads TO {num_threads}")
    con.execute("SET preserve_insertion_order = false")  # Allow sort optimization

    # Configure memory and temp directory for large files
    if temp_dir:
        con.execute(f"SET temp_directory = '{temp_dir}'")

    # DuckDB automatically uses external sorting when needed
    # No explicit setting required in recent versions

    start_time = time.time()

    # Get column names from input file
    print("Reading schema...")
    schema_result = con.execute(f"DESCRIBE SELECT * FROM '{input_file}'").fetchall()
    columns = [row[0] for row in schema_result]
    print(f"Columns: {len(columns)}")

    # Check if we have ref_id and position columns
    has_ref_id = 'ref_id' in columns
    has_position = 'position' in columns

    if not has_position:
        raise ValueError("Input file must have 'position' column for sorting")

    # Build sort key
    if has_ref_id:
        sort_key = "ref_id, position"
        print(f"Sorting by: ref_id, position")
    else:
        sort_key = "position"
        print(f"Sorting by: position (no ref_id column)")

    # Count records
    print("Counting records...")
    count_start = time.time()
    record_count = con.execute(f"SELECT COUNT(*) FROM '{input_file}'").fetchone()[0]
    count_time = time.time() - count_start
    print(f"Records: {record_count:,} (counted in {count_time:.1f}s)")

    # Build COPY query with sort
    # DuckDB will use external merge sort for large files
    print()
    print("Sorting and writing...")
    sort_start = time.time()

    # Use COPY with ORDER BY for sorted output
    query = f"""
        COPY (
            SELECT *
            FROM '{input_file}'
            ORDER BY {sort_key}
        )
        TO '{output_file}'
        (
            FORMAT PARQUET,
            COMPRESSION 'ZSTD',
            COMPRESSION_LEVEL {compression_level},
            ROW_GROUP_SIZE {row_group_size}
        )
    """

    con.execute(query)
    sort_time = time.time() - sort_start
    total_time = time.time() - start_time

    # Get output size
    output_size = os.path.getsize(output_file)

    # Calculate statistics
    size_diff = input_size - output_size
    size_pct = (size_diff / input_size * 100) if input_size > 0 else 0

    print()
    print("=" * 70)
    print("SORTING COMPLETE")
    print("=" * 70)
    print(f"Sort+write time: {sort_time:.1f}s")
    print(f"Total time: {total_time:.1f}s")
    print(f"Throughput: {record_count / sort_time:,.0f} records/sec")
    print()
    print(f"Input size:  {input_size / 1024 / 1024:.1f} MB")
    print(f"Output size: {output_size / 1024 / 1024:.1f} MB")
    print(f"Difference:  {size_diff / 1024 / 1024:+.1f} MB ({size_pct:+.1f}%)")
    print()

    con.close()

    return {
        'input_file': input_file,
        'output_file': output_file,
        'record_count': record_count,
        'input_size_bytes': input_size,
        'output_size_bytes': output_size,
        'size_reduction_bytes': size_diff,
        'size_reduction_pct': size_pct,
        'sort_time_seconds': sort_time,
        'total_time_seconds': total_time,
        'throughput_records_per_sec': record_count / sort_time if sort_time > 0 else 0,
    }


def compare_compression_detail(unsorted_file: str, sorted_file: str):
    """Compare per-column compression between unsorted and sorted files."""
    import pyarrow.parquet as pq

    def get_column_sizes(filename):
        pf = pq.ParquetFile(filename)
        meta = pf.metadata
        sizes = {}
        for i in range(meta.num_row_groups):
            rg = meta.row_group(i)
            for j in range(rg.num_columns):
                col = rg.column(j)
                name = col.path_in_schema
                if name not in sizes:
                    sizes[name] = {'compressed': 0, 'uncompressed': 0}
                sizes[name]['compressed'] += col.total_compressed_size
                sizes[name]['uncompressed'] += col.total_uncompressed_size
        return sizes

    print("\nPer-column compression comparison:")
    print(f"{'Column':<25} {'Unsorted':>12} {'Sorted':>12} {'Savings':>12} {'%':>8}")
    print("-" * 72)

    unsorted_sizes = get_column_sizes(unsorted_file)
    sorted_sizes = get_column_sizes(sorted_file)

    # Key columns to highlight
    key_columns = ['ref_id', 'position', 'aligned_length', 'is_reverse', 'mapq',
                   'MD', 'quality', 'sequence_packed', 'cigar', 'AS', 'NM']

    for col in key_columns:
        if col in unsorted_sizes and col in sorted_sizes:
            u = unsorted_sizes[col]['compressed']
            s = sorted_sizes[col]['compressed']
            savings = u - s
            pct = (savings / u * 100) if u > 0 else 0
            print(f"{col:<25} {u/1024/1024:>10.2f}MB {s/1024/1024:>10.2f}MB {savings/1024/1024:>+10.2f}MB {pct:>+7.1f}%")

    # Total
    total_u = sum(s['compressed'] for s in unsorted_sizes.values())
    total_s = sum(s['compressed'] for s in sorted_sizes.values())
    savings = total_u - total_s
    pct = (savings / total_u * 100) if total_u > 0 else 0
    print("-" * 72)
    print(f"{'TOTAL':<25} {total_u/1024/1024:>10.2f}MB {total_s/1024/1024:>10.2f}MB {savings/1024/1024:>+10.2f}MB {pct:>+7.1f}%")


def verify_sort_order(parquet_file: str, sample_size: int = 1000) -> bool:
    """Verify that a Parquet file is sorted by (ref_id, position)."""
    import pyarrow.parquet as pq

    print(f"\nVerifying sort order (first {sample_size} rows)...")

    # Read sample
    table = pq.read_table(parquet_file, columns=['ref_id', 'position'])

    if table.num_rows == 0:
        print("  Empty file - OK")
        return True

    # Get values
    ref_ids = table['ref_id'].to_pylist()[:sample_size]
    positions = table['position'].to_pylist()[:sample_size]

    # Check order
    for i in range(1, len(ref_ids)):
        if (ref_ids[i], positions[i]) < (ref_ids[i-1], positions[i-1]):
            print(f"  VIOLATION at row {i}:")
            print(f"    ({ref_ids[i-1]}, {positions[i-1]}) > ({ref_ids[i]}, {positions[i]})")
            return False

    print(f"  Sort order verified OK")
    return True


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Sort Parquet file by (ref_id, position) using DuckDB"
    )
    parser.add_argument('input_file', help="Input Parquet file")
    parser.add_argument('output_file', help="Output sorted Parquet file")
    parser.add_argument('--compression-level', '-c', type=int, default=6,
                        help="ZSTD compression level (default: 6)")
    parser.add_argument('--row-group-size', '-r', type=int, default=500000,
                        help="Row group size (default: 500000)")
    parser.add_argument('--threads', '-t', type=int, default=8,
                        help="Number of threads (default: 8)")
    parser.add_argument('--temp-dir', help="Temporary directory for disk spilling")
    parser.add_argument('--compare', action='store_true',
                        help="Compare per-column compression")
    parser.add_argument('--verify', action='store_true',
                        help="Verify sort order after sorting")

    args = parser.parse_args()

    # Run sorting
    stats = sort_parquet_by_position(
        input_file=args.input_file,
        output_file=args.output_file,
        compression_level=args.compression_level,
        row_group_size=args.row_group_size,
        num_threads=args.threads,
        temp_dir=args.temp_dir,
    )

    # Optional: compare compression
    if args.compare:
        compare_compression_detail(args.input_file, args.output_file)

    # Optional: verify sort
    if args.verify:
        verify_sort_order(args.output_file)

    print("\nDone!")


if __name__ == "__main__":
    main()
