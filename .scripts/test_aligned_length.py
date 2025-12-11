#!/usr/bin/env python3
"""
Test aligned_length column compression vs end_position.

Based on model consensus:
- aligned_length (uint16) clusters around 100-150bp, compresses much better
- end_position (int32) has high entropy like position, compresses poorly
"""

import sys
sys.path.insert(0, '/maps/projects/fernandezguerra/apps/repos/bam-filter')

from bam_filter.parquet_converter_optimized import convert_bam_to_optimized_parquet
import pyarrow.parquet as pq
import os


def test_aligned_length():
    """Test aligned_length column compression."""
    input_file = "/scratch/tmp/zs_test/test_small.sam.gz"
    output_file = "/scratch/tmp/zs_test/test_aligned_length.parquet"

    # Remove existing output
    if os.path.exists(output_file):
        os.remove(output_file)

    print("Testing aligned_length column (replaces end_position)...")
    stats = convert_bam_to_optimized_parquet(
        input_file=input_file,
        output_file=output_file,
        batch_size=10000,
        compression_level=3,
        num_threads=1,
        calculate_pmd=False,
        library_type="ds",
        store_sequences=True
    )

    print(f"\nProcessed {stats['total_records']} records in {stats['processing_time_seconds']:.2f}s")

    # Check schema
    table = pq.read_table(output_file)
    print(f"\nSchema columns: {table.schema.names}")

    # Verify columns
    assert 'aligned_length' in table.schema.names, "Missing aligned_length column"
    assert 'is_reverse' in table.schema.names, "Missing is_reverse column"
    assert 'end_position' not in table.schema.names, "end_position should be removed"
    print("Schema correct: aligned_length present, end_position removed.")

    # Get column sizes
    pf = pq.ParquetFile(output_file)
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

    # Compare relevant columns
    print(f"\nColumn sizes:")
    print(f"{'Column':<20} {'Compressed':>12} {'Uncompressed':>14} {'Ratio':>8}")
    print("-" * 58)
    for col in ['position', 'aligned_length', 'is_reverse']:
        s = sizes[col]
        ratio = s['compressed'] / s['uncompressed'] if s['uncompressed'] > 0 else 0
        print(f"{col:<20} {s['compressed']/1024/1024:>10.2f}MB {s['uncompressed']/1024/1024:>12.2f}MB {ratio:>7.1%}")

    # File size comparison
    file_size = os.path.getsize(output_file)
    nrows = table.num_rows
    bytes_per_rec = file_size / nrows

    print(f"\nFile stats:")
    print(f"  Total size: {file_size/1024/1024:.1f} MB")
    print(f"  Records: {nrows:,}")
    print(f"  Bytes/record: {bytes_per_rec:.2f}")

    # Sample data
    aligned_lengths = table['aligned_length'].to_pylist()[:1000]
    is_reverse = table['is_reverse'].to_pylist()[:1000]
    positions = table['position'].to_pylist()[:1000]

    print(f"\nAligned length distribution (first 1000):")
    from collections import Counter
    length_counts = Counter(aligned_lengths)
    top_lengths = length_counts.most_common(10)
    print(f"  Top 10 lengths: {top_lengths}")
    print(f"  Min: {min(aligned_lengths)}, Max: {max(aligned_lengths)}, Mean: {sum(aligned_lengths)/len(aligned_lengths):.1f}")

    # Verify end_position can be computed
    print(f"\nSample (pos + aligned_length = end_position):")
    for i in range(5):
        end = positions[i] + aligned_lengths[i]
        print(f"  {positions[i]} + {aligned_lengths[i]} = {end}")

    print("\nTest completed successfully!")


if __name__ == "__main__":
    test_aligned_length()
