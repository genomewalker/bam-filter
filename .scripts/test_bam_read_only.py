#!/usr/bin/env python3
"""Benchmark BAM reading speed using existing Cython HTSlib bindings."""

import sys
import time
from bam_filter.parquet_converter_pure_cpp import convert_sam_bam_to_parquet_pure_cpp

def main():
    if len(sys.argv) < 2:
        print("Usage: python test_bam_read_only.py <input.bam>")
        sys.exit(1)

    input_file = sys.argv[1]

    # First test: no processing at all, just iterate
    print("=" * 60)
    print("Testing conversion with MINIMAL schema")
    print("=" * 60)
    print(f"Input: {input_file}")
    print()

    # We'll test with different batch sizes to see if that's the bottleneck
    for batch_size in [100000, 200000, 500000]:
        print(f"\n--- Batch size: {batch_size:,} ---")
        import tempfile
        import shutil

        output_dir = tempfile.mkdtemp(prefix="bam_bench_")

        start_time = time.time()

        stats = convert_sam_bam_to_parquet_pure_cpp(
            input_file,
            output_dir,
            num_partitions=16,
            batch_size=batch_size,
            write_by_reference=True,
            write_by_read=False,
            compression_level=1,  # Fastest compression
            num_threads=8,
            calculate_pmd=False,
            library_type="ds",
        )

        elapsed = time.time() - start_time
        records = stats["total_records"]
        rate = records / elapsed if elapsed > 0 else 0

        print(f"Records: {records:,}")
        print(f"Time: {elapsed:.1f}s")
        print(f"Throughput: {rate:,.0f} rec/s ({rate/1e6:.3f} M/s)")

        shutil.rmtree(output_dir)

if __name__ == "__main__":
    main()
