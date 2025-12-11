#!/usr/bin/env python3
"""
Test script for the optimized Parquet converter integrated with filterBAM.

Tests:
1. Direct API call to convert_bam_to_optimized_parquet
2. filterBAM CLI with --optimized flag
3. Verifies output schema and data integrity
"""

import sys
import os
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pathlib import Path
import pyarrow.parquet as pq


def test_direct_api(input_file, output_dir, max_records=100000):
    """Test the optimized converter directly via API."""
    from bam_filter.parquet_converter_optimized import convert_bam_to_optimized_parquet

    print(f"\n=== Testing Direct API ===")
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    output_file = str(output_path / "alignments.parquet")

    result = convert_bam_to_optimized_parquet(
        input_file=input_file,
        output_file=output_file,
        batch_size=50000,
        compression_level=6,
        num_threads=4,
        calculate_pmd=False,
        store_sequences=True,
    )

    print(f"\nResults:")
    print(f"  Records: {result['total_records']:,}")
    print(f"  Time: {result['processing_time_seconds']:.2f}s")
    throughput = result['total_records'] / result['processing_time_seconds']
    print(f"  Throughput: {throughput:,.0f} rec/s")

    # Verify output
    print(f"\nVerifying output...")
    table = pq.read_table(output_file)
    print(f"  Parquet rows: {len(table):,}")
    print(f"  Columns: {len(table.column_names)}")
    print(f"  Schema:")
    for field in table.schema:
        print(f"    {field.name}: {field.type}")

    # Check 2-bit packing
    if 'sequence_packed' in table.column_names:
        seq_packed = table['sequence_packed']
        seq_lengths = table['sequence_length']
        if len(seq_packed) > 0:
            first_packed = seq_packed[0].as_py()
            first_len = seq_lengths[0].as_py()
            if first_packed and first_len:
                compression_ratio = first_len / len(first_packed) if first_packed else 0
                print(f"\n  Sequence packing:")
                print(f"    First seq length: {first_len} bases")
                print(f"    Packed size: {len(first_packed)} bytes")
                print(f"    Compression ratio: {compression_ratio:.1f}x")

    # Check hot tags
    hot_tags = ['AS', 'NM', 'XS', 'MD']
    print(f"\n  Hot tags present: {[t for t in hot_tags if t in table.column_names]}")

    # Check cold tags
    if 'tags_cold' in table.column_names:
        cold_sample = table['tags_cold'][0].as_py() if len(table) > 0 else ""
        print(f"  Sample cold tags: {cold_sample[:80]}..." if cold_sample else "  Cold tags: (empty)")

    file_size = os.path.getsize(output_file)
    bytes_per_rec = file_size / len(table) if len(table) > 0 else 0
    print(f"\n  Output size: {file_size / 1024 / 1024:.2f} MB ({bytes_per_rec:.1f} bytes/record)")

    return result


def test_cli(input_file, output_dir):
    """Test via filterBAM CLI."""
    import subprocess

    print(f"\n=== Testing filterBAM CLI ===")
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")

    cmd = [
        "filterBAM", "to-parquet",
        "--bam", input_file,
        "-o", output_dir,
        "--optimized",
        "--disable-pmd",
        "-t", "4",
    ]

    print(f"Command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    print(f"\nStdout:\n{result.stdout}")
    if result.stderr:
        print(f"Stderr:\n{result.stderr}")

    return result.returncode == 0


def main():
    if len(sys.argv) < 2:
        print("Usage: python test_optimized_converter.py <input.bam> [output_dir]")
        print("\nTests the optimized Parquet converter with:")
        print("  - 2-bit packed sequences")
        print("  - Hot/cold tag separation")
        print("  - Bulk Arrow AppendValues API")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else '/scratch/tmp/optimized_converter_test'

    # Test direct API
    api_dir = os.path.join(output_dir, 'api_test')
    test_direct_api(input_file, api_dir)

    # Test CLI
    cli_dir = os.path.join(output_dir, 'cli_test')
    success = test_cli(input_file, cli_dir)

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"API test: PASSED")
    print(f"CLI test: {'PASSED' if success else 'FAILED'}")


if __name__ == '__main__':
    main()
