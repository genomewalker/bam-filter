#!/usr/bin/env python3
"""
Convert BAM/SAM to sorted Parquet with optimal compression.

Two-phase approach:
1. Convert BAM → Parquet (streaming, uses C++ Arrow writer)
2. Sort by read_id with DuckDB (external merge sort, memory efficient)

Sorting by read_id provides:
- 39% better compression (multi-mapped reads cluster together)
- Optimal for LCA and EM algorithms (group by read_id)
- Preserves natural read property clustering within each read
"""

import os
import sys
import time
import tempfile
from pathlib import Path
from typing import Dict, Any, Optional

sys.path.insert(0, '/maps/projects/fernandezguerra/apps/repos/bam-filter')


def convert_and_sort(
    input_file: str,
    output_file: str,
    batch_size: int = 10000,
    compression_level: int = 6,
    num_threads: int = 4,
    calculate_pmd: bool = False,
    library_type: str = "ds",
    store_sequences: bool = True,
    row_group_size: int = 100000,
    keep_unsorted: bool = False,
) -> Dict[str, Any]:
    """
    Convert BAM/SAM to sorted Parquet with optimal compression.

    Args:
        input_file: Input BAM/SAM file
        output_file: Output Parquet file (sorted by read_id)
        batch_size: Batch size for conversion
        compression_level: ZSTD compression level (1-22)
        num_threads: Number of threads for conversion and sorting
        calculate_pmd: Calculate PMD scores during conversion
        library_type: Library type for PMD calculation ("ss" or "ds")
        store_sequences: Store sequences and quality scores
        row_group_size: Parquet row group size
        keep_unsorted: Keep unsorted intermediate file

    Returns:
        Dictionary with conversion statistics
    """
    import duckdb
    from bam_filter.parquet_converter_optimized import convert_bam_to_optimized_parquet

    print("=" * 70)
    print("BAM/SAM to Sorted Parquet Converter")
    print("=" * 70)
    print(f"Input:  {input_file}")
    print(f"Output: {output_file}")
    print(f"Threads: {num_threads}")
    print(f"Compression level: {compression_level}")
    print(f"Row group size: {row_group_size:,}")
    print()

    # Phase 1: Convert to unsorted Parquet
    output_dir = Path(output_file).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    # Use temp file for unsorted output
    unsorted_file = str(output_file) + ".unsorted"

    print("Phase 1: Converting BAM/SAM to Parquet...")
    convert_start = time.time()

    convert_stats = convert_bam_to_optimized_parquet(
        input_file=input_file,
        output_file=unsorted_file,
        batch_size=batch_size,
        compression_level=compression_level,
        num_threads=num_threads,
        calculate_pmd=calculate_pmd,
        library_type=library_type,
        store_sequences=store_sequences,
    )

    convert_time = time.time() - convert_start
    unsorted_size = os.path.getsize(unsorted_file)

    print(f"  Records: {convert_stats['total_records']:,}")
    print(f"  Time: {convert_time:.1f}s")
    print(f"  Throughput: {convert_stats['total_records']/convert_time:,.0f} rec/s")
    print(f"  Unsorted size: {unsorted_size/1e6:.1f} MB")
    print()

    # Phase 2: Sort by read_id with DuckDB
    print("Phase 2: Sorting by read_id with DuckDB...")
    sort_start = time.time()

    con = duckdb.connect(":memory:")
    con.execute(f"SET threads TO {num_threads}")

    query = f"""
        COPY (
            SELECT *
            FROM read_parquet('{unsorted_file}')
            ORDER BY read_id
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
    con.close()

    sort_time = time.time() - sort_start
    sorted_size = os.path.getsize(output_file)

    print(f"  Time: {sort_time:.1f}s")
    print(f"  Throughput: {convert_stats['total_records']/sort_time:,.0f} rec/s")
    print(f"  Sorted size: {sorted_size/1e6:.1f} MB")
    print()

    # Cleanup unsorted file
    if not keep_unsorted:
        os.remove(unsorted_file)
        print(f"  Removed unsorted file")
    else:
        print(f"  Kept unsorted file: {unsorted_file}")

    # Summary
    total_time = convert_time + sort_time
    size_savings = unsorted_size - sorted_size
    size_pct = (size_savings / unsorted_size * 100) if unsorted_size > 0 else 0

    print()
    print("=" * 70)
    print("CONVERSION COMPLETE")
    print("=" * 70)
    print(f"Total records: {convert_stats['total_records']:,}")
    print(f"Total time: {total_time:.1f}s ({total_time/60:.1f} min)")
    print(f"  Convert: {convert_time:.1f}s")
    print(f"  Sort:    {sort_time:.1f}s")
    print(f"Overall throughput: {convert_stats['total_records']/total_time:,.0f} rec/s")
    print()
    print(f"Size comparison:")
    print(f"  Unsorted: {unsorted_size/1e6:.1f} MB")
    print(f"  Sorted:   {sorted_size/1e6:.1f} MB")
    print(f"  Savings:  {size_savings/1e6:.1f} MB ({size_pct:.1f}%)")
    print()

    return {
        'total_records': convert_stats['total_records'],
        'convert_time_seconds': convert_time,
        'sort_time_seconds': sort_time,
        'total_time_seconds': total_time,
        'unsorted_size_bytes': unsorted_size,
        'sorted_size_bytes': sorted_size,
        'size_savings_bytes': size_savings,
        'size_savings_pct': size_pct,
        'throughput_overall': convert_stats['total_records'] / total_time,
    }


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert BAM/SAM to sorted Parquet with optimal compression"
    )
    parser.add_argument('input_file', help="Input BAM/SAM file")
    parser.add_argument('output_file', help="Output Parquet file")
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help="Number of threads (default: 4)")
    parser.add_argument('-c', '--compression-level', type=int, default=6,
                        help="ZSTD compression level (default: 6)")
    parser.add_argument('-b', '--batch-size', type=int, default=10000,
                        help="Batch size (default: 10000)")
    parser.add_argument('-r', '--row-group-size', type=int, default=100000,
                        help="Row group size (default: 100000)")
    parser.add_argument('--pmd', action='store_true',
                        help="Calculate PMD scores")
    parser.add_argument('--library-type', choices=['ss', 'ds'], default='ds',
                        help="Library type for PMD (default: ds)")
    parser.add_argument('--no-sequences', action='store_true',
                        help="Don't store sequences and quality scores")
    parser.add_argument('--keep-unsorted', action='store_true',
                        help="Keep unsorted intermediate file")

    args = parser.parse_args()

    stats = convert_and_sort(
        input_file=args.input_file,
        output_file=args.output_file,
        batch_size=args.batch_size,
        compression_level=args.compression_level,
        num_threads=args.threads,
        calculate_pmd=args.pmd,
        library_type=args.library_type,
        store_sequences=not args.no_sequences,
        row_group_size=args.row_group_size,
        keep_unsorted=args.keep_unsorted,
    )

    print("Done!")


if __name__ == "__main__":
    main()
