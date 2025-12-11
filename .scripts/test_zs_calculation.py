#!/usr/bin/env python3
"""
Test ZS score calculation using Cython functions directly.
"""

import sys
sys.path.insert(0, '/maps/projects/fernandezguerra/apps/repos/bam-filter')

from bam_filter.parquet_converter_optimized import convert_bam_to_optimized_parquet
import pyarrow.parquet as pq


def test_zs_single_file():
    """Test ZS calculation on a small file."""
    input_file = "/scratch/tmp/zs_test/test_small.sam.gz"
    output_file = "/scratch/tmp/zs_test/test_zs_single.parquet"

    print("Testing single-file optimized converter...")
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

    print(f"Processed {stats['total_records']} records in {stats['processing_time_seconds']:.2f}s")

    # Check output
    table = pq.read_table(output_file)
    zs_scores = table['zs_score'].to_pylist()

    non_null = [x for x in zs_scores if x > -1e19]
    print(f"\nZS Score Results:")
    print(f"  Total records: {len(zs_scores)}")
    print(f"  Records with ZS: {len(non_null)}")
    if non_null:
        print(f"  Sample ZS: {non_null[:5]}")
        print(f"  Min ZS: {min(non_null):.4f}")
        print(f"  Max ZS: {max(non_null):.4f}")

        # Check for non-zero values
        non_zero = [x for x in non_null if abs(x) > 1e-6]
        print(f"  Non-zero ZS: {len(non_zero)}")
        if non_zero:
            print(f"  Sample non-zero: {non_zero[:5]}")


if __name__ == "__main__":
    test_zs_single_file()
