#!/usr/bin/env python3
"""
Test script for BAM to Parquet conversion using the Cython/C++ implementation.

This script demonstrates and tests the high-performance BAM to Parquet converter
that uses HTSlib for reading and Arrow C++ for writing.

Usage:
    python .scripts/test_bam_to_parquet.py <input.bam> <output_dir> [--threads N]

Example:
    python .scripts/test_bam_to_parquet.py sample.bam /tmp/parquet_test --threads 8
"""

import argparse
import sys
import time
from pathlib import Path


def test_single_worker(input_file: str, output_dir: str, batch_size: int = 100000):
    """Test single-worker conversion (useful for debugging)."""
    from bam_filter.parquet_converter_pure_cpp import (
        convert_sam_bam_to_parquet_pure_cpp,
        save_bam_header,
    )

    print(f"\n{'='*60}")
    print("SINGLE WORKER TEST")
    print(f"{'='*60}")
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print(f"Batch size: {batch_size:,}")
    print()

    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Save BAM header
    print("Saving BAM header...")
    save_bam_header(input_file, output_dir)

    # Convert
    print("\nStarting conversion...")
    start_time = time.time()

    stats = convert_sam_bam_to_parquet_pure_cpp(
        input_file,
        output_dir,
        num_partitions=16,
        batch_size=batch_size,
        write_by_reference=True,
        write_by_read=True,
        compression_level=3,
        num_threads=4,
        skip_records=0,
        max_records=0,  # 0 = all records
        calculate_pmd=True,
        library_type="ds",
    )

    elapsed = time.time() - start_time
    records = stats["total_records"]
    rate = records / elapsed if elapsed > 0 else 0

    print(f"\n{'='*60}")
    print("RESULTS")
    print(f"{'='*60}")
    print(f"Records processed: {records:,}")
    print(f"Time: {elapsed:.2f}s")
    print(f"Throughput: {rate:,.0f} records/sec ({rate/1e6:.3f} M/sec)")

    # List output files
    output_path = Path(output_dir)
    parquet_files = list(output_path.rglob("*.parquet"))
    total_size = sum(f.stat().st_size for f in parquet_files)
    print(f"\nOutput files: {len(parquet_files)}")
    print(f"Total size: {total_size / 1e6:.2f} MB")

    return stats


def test_parallel_workers(
    input_file: str, output_dir: str, num_processes: int = 8, batch_size: int = 100000
):
    """Test parallel conversion using ProcessPoolExecutor."""
    from bam_filter.parallel_parquet_records import parallel_convert_records

    print(f"\n{'='*60}")
    print("PARALLEL WORKERS TEST")
    print(f"{'='*60}")
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print(f"Processes: {num_processes}")
    print(f"Batch size: {batch_size:,}")
    print()

    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    start_time = time.time()

    stats = parallel_convert_records(
        input_file=input_file,
        output_dir=output_dir,
        num_processes=num_processes,
        num_partitions=16,
        batch_size=batch_size,
        compression_level=3,
        num_threads_per_worker=1,
        count_records=False,
        calculate_pmd=True,
        library_type="ds",
    )

    elapsed = time.time() - start_time
    records = stats["total_records"]
    rate = records / elapsed if elapsed > 0 else 0

    print(f"\n{'='*60}")
    print("PARALLEL RESULTS")
    print(f"{'='*60}")
    print(f"Records processed: {records:,}")
    print(f"Time: {elapsed:.2f}s")
    print(f"Throughput: {rate:,.0f} records/sec ({rate/1e6:.3f} M/sec)")

    # List output files
    output_path = Path(output_dir)
    parquet_files = list(output_path.rglob("*.parquet"))
    total_size = sum(f.stat().st_size for f in parquet_files)
    print(f"\nOutput files: {len(parquet_files)}")
    print(f"Total size: {total_size / 1e6:.2f} MB")

    return stats


def verify_parquet_output(output_dir: str, sample_rows: int = 5):
    """Verify Parquet output by reading sample rows."""
    try:
        import pyarrow.parquet as pq
    except ImportError:
        print("pyarrow not available for verification")
        return

    print(f"\n{'='*60}")
    print("VERIFICATION")
    print(f"{'='*60}")

    output_path = Path(output_dir)

    # Find first parquet file
    parquet_files = list(output_path.rglob("*.parquet"))
    if not parquet_files:
        print("No Parquet files found!")
        return

    sample_file = parquet_files[0]
    print(f"Sample file: {sample_file}")

    # Read metadata
    pf = pq.ParquetFile(sample_file)
    print(f"Schema: {pf.schema_arrow}")
    print(f"Num row groups: {pf.metadata.num_row_groups}")
    print(f"Num rows: {pf.metadata.num_rows}")

    # Read sample
    table = pf.read()
    df = table.to_pandas()
    print(f"\nSample rows ({sample_rows}):")
    print(df.head(sample_rows).to_string())


def main():
    parser = argparse.ArgumentParser(
        description="Test BAM to Parquet conversion",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("input_file", help="Input BAM/SAM file")
    parser.add_argument("output_dir", help="Output directory for Parquet files")
    parser.add_argument(
        "--threads", type=int, default=8, help="Number of parallel processes (default: 8)"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=100000,
        help="Batch size for writing (default: 100000)",
    )
    parser.add_argument(
        "--single",
        action="store_true",
        help="Run single-worker test only",
    )
    parser.add_argument(
        "--parallel",
        action="store_true",
        help="Run parallel test only",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Verify output by reading sample rows",
    )
    parser.add_argument(
        "--max-records",
        type=int,
        default=0,
        help="Maximum records to process (0 = all, useful for quick tests)",
    )

    args = parser.parse_args()

    # Validate input
    if not Path(args.input_file).exists():
        print(f"Error: Input file not found: {args.input_file}")
        sys.exit(1)

    run_single = args.single or (not args.single and not args.parallel)
    run_parallel = args.parallel or (not args.single and not args.parallel)

    # Run tests
    if run_single:
        single_output = f"{args.output_dir}/single_worker"
        test_single_worker(args.input_file, single_output, args.batch_size)
        if args.verify:
            verify_parquet_output(single_output)

    if run_parallel and not args.max_records:
        parallel_output = f"{args.output_dir}/parallel"
        test_parallel_workers(
            args.input_file, parallel_output, args.threads, args.batch_size
        )
        if args.verify:
            verify_parquet_output(parallel_output)

    print("\nDone!")


if __name__ == "__main__":
    main()
