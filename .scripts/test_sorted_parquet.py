#!/usr/bin/env python3
"""
Test compression improvement from sorting Parquet by (ref_id, position).

Sorting clusters similar alignments together, improving:
1. RLE compression for ref_id (long runs per chromosome)
2. Delta compression for position (smaller deltas)
3. ZSTD compression for MD/quality (similar patterns adjacent)
"""

import sys
sys.path.insert(0, '/maps/projects/fernandezguerra/apps/repos/bam-filter')

import pyarrow.parquet as pq
import pyarrow.compute as pc
import pyarrow as pa
import os
import time


def sort_parquet_by_position(input_file: str, output_file: str, compression_level: int = 6):
    """
    Sort a Parquet file by (ref_id, position) for better compression.

    Args:
        input_file: Input Parquet file
        output_file: Output sorted Parquet file
        compression_level: ZSTD compression level
    """
    print(f"Reading {input_file}...")
    start = time.time()
    table = pq.read_table(input_file)
    read_time = time.time() - start
    print(f"  Read {table.num_rows:,} rows in {read_time:.1f}s")

    print(f"Sorting by (ref_id, position)...")
    start = time.time()
    # Sort by ref_id, then position
    sort_indices = pc.sort_indices(table, sort_keys=[("ref_id", "ascending"), ("position", "ascending")])
    sorted_table = table.take(sort_indices)
    sort_time = time.time() - start
    print(f"  Sorted in {sort_time:.1f}s")

    print(f"Writing sorted Parquet to {output_file}...")
    start = time.time()
    pq.write_table(
        sorted_table,
        output_file,
        compression='zstd',
        compression_level=compression_level,
        use_dictionary=['ref_id', 'cigar', 'read_group'],
        row_group_size=500000,  # ~500K rows per row group
    )
    write_time = time.time() - start
    print(f"  Wrote in {write_time:.1f}s")

    return {
        'read_time': read_time,
        'sort_time': sort_time,
        'write_time': write_time,
    }


def compare_compression(unsorted_file: str, sorted_file: str):
    """Compare column sizes between unsorted and sorted files."""

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

    unsorted_sizes = get_column_sizes(unsorted_file)
    sorted_sizes = get_column_sizes(sorted_file)

    print("\nColumn size comparison (unsorted vs sorted):")
    print(f"{'Column':<20} {'Unsorted':>12} {'Sorted':>12} {'Savings':>10} {'%':>8}")
    print("-" * 65)

    total_unsorted = 0
    total_sorted = 0

    # Key columns to highlight
    key_columns = ['ref_id', 'position', 'MD', 'quality', 'sequence_packed', 'aligned_length', 'is_reverse']

    for col in key_columns:
        if col in unsorted_sizes and col in sorted_sizes:
            u = unsorted_sizes[col]['compressed']
            s = sorted_sizes[col]['compressed']
            savings = u - s
            pct = (savings / u * 100) if u > 0 else 0
            total_unsorted += u
            total_sorted += s
            print(f"{col:<20} {u/1024/1024:>10.2f}MB {s/1024/1024:>10.2f}MB {savings/1024/1024:>+9.2f}MB {pct:>+7.1f}%")

    # Total for all columns
    total_u = sum(s['compressed'] for s in unsorted_sizes.values())
    total_s = sum(s['compressed'] for s in sorted_sizes.values())
    savings = total_u - total_s
    pct = (savings / total_u * 100) if total_u > 0 else 0

    print("-" * 65)
    print(f"{'TOTAL':<20} {total_u/1024/1024:>10.2f}MB {total_s/1024/1024:>10.2f}MB {savings/1024/1024:>+9.2f}MB {pct:>+7.1f}%")

    # File sizes
    unsorted_size = os.path.getsize(unsorted_file)
    sorted_size = os.path.getsize(sorted_file)

    print(f"\nFile sizes:")
    print(f"  Unsorted: {unsorted_size/1024/1024:.1f} MB")
    print(f"  Sorted:   {sorted_size/1024/1024:.1f} MB")
    print(f"  Savings:  {(unsorted_size-sorted_size)/1024/1024:.1f} MB ({(unsorted_size-sorted_size)/unsorted_size*100:.1f}%)")


def main():
    unsorted_file = "/scratch/tmp/zs_test/test_aligned_length.parquet"
    sorted_file = "/scratch/tmp/zs_test/test_sorted.parquet"

    if not os.path.exists(unsorted_file):
        print(f"Error: {unsorted_file} not found")
        return

    # Remove existing output
    if os.path.exists(sorted_file):
        os.remove(sorted_file)

    print("=" * 65)
    print("Testing compression improvement from position sorting")
    print("=" * 65)

    # Sort and write
    stats = sort_parquet_by_position(unsorted_file, sorted_file, compression_level=6)

    print(f"\nTiming summary:")
    print(f"  Read:  {stats['read_time']:.1f}s")
    print(f"  Sort:  {stats['sort_time']:.1f}s")
    print(f"  Write: {stats['write_time']:.1f}s")
    print(f"  Total: {sum(stats.values()):.1f}s")

    # Compare
    compare_compression(unsorted_file, sorted_file)

    # Verify sorted order
    print("\nVerifying sort order...")
    sorted_table = pq.read_table(sorted_file, columns=['ref_id', 'position'])
    ref_ids = sorted_table['ref_id'].to_pylist()[:100]
    positions = sorted_table['position'].to_pylist()[:100]

    is_sorted = True
    for i in range(1, len(ref_ids)):
        if (ref_ids[i], positions[i]) < (ref_ids[i-1], positions[i-1]):
            is_sorted = False
            print(f"  Sort violation at {i}: ({ref_ids[i-1]}, {positions[i-1]}) > ({ref_ids[i]}, {positions[i]})")
            break

    if is_sorted:
        print("  Sort order verified OK")

    print("\nDone!")


if __name__ == "__main__":
    main()
