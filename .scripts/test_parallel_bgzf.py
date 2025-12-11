#!/usr/bin/env python3
"""
Test parallel BGZF block-based BAM/SAM.gz to Parquet conversion.

This tests the optimized parallel converter that:
1. Scans BGZF block boundaries
2. Spawns parallel workers that seek directly to byte ranges
3. Each worker uses the optimized Cython extract_record_nogil()
4. Writes to separate Parquet files (can be merged with DuckDB/polars)

Usage:
    python test_parallel_bgzf.py <input.bam|sam.gz> <output_dir> [num_workers]
"""

import sys
import time
from pathlib import Path


def main():
    if len(sys.argv) < 3:
        print("Usage: python test_parallel_bgzf.py <input.bam|sam.gz> <output_dir> [num_workers]")
        print("\nExamples:")
        print("  python test_parallel_bgzf.py input.bam /tmp/output 8")
        print("  python test_parallel_bgzf.py sample.sam.gz /tmp/output 4")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    num_workers = int(sys.argv[3]) if len(sys.argv) > 3 else 8

    # Import the parallel converter
    from bam_filter.parallel_bgzf import parallel_convert_bgzf

    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print(f"Workers: {num_workers}")
    print(f"Options: NO PMD, NO SEQUENCES (max speed)")
    print()

    # Clean output directory
    output_path = Path(output_dir)
    if output_path.exists():
        import shutil
        shutil.rmtree(output_path)

    start_time = time.time()

    stats = parallel_convert_bgzf(
        input_file,
        output_dir,
        num_workers=num_workers,
        batch_size=100000,
        compression_level=3,
        num_threads_per_worker=1,
        calculate_pmd=False,
        library_type="ds",
        store_sequences=False,
    )

    total_elapsed = time.time() - start_time
    total_records = stats['total_records']

    print(f"\n{'='*60}")
    print(f"FINAL RESULTS")
    print(f"{'='*60}")
    print(f"Total records: {total_records:,}")
    print(f"Total time: {total_elapsed:.1f}s")
    print(f"Throughput: {total_records/total_elapsed:,.0f} rec/s")
    print(f"           {total_records/total_elapsed/1e6:.3f} M/s")
    print()

    # Output size
    parquet_files = list(Path(output_dir).rglob("*.parquet"))
    total_size = sum(f.stat().st_size for f in parquet_files)
    print(f"Output: {len(parquet_files)} Parquet files, {total_size/1e6:.1f} MB")

    # Show per-worker stats
    print(f"\nPer-worker breakdown:")
    for r in sorted(stats['worker_results'], key=lambda x: x['worker_id']):
        worker_id = r['worker_id']
        records = r['total_records']
        elapsed = r['elapsed']
        rate = records / elapsed if elapsed > 0 else 0
        if 'error' in r:
            print(f"  Worker {worker_id}: ERROR - {r['error']}")
        else:
            print(f"  Worker {worker_id}: {records:>10,} records in {elapsed:>5.1f}s ({rate:>8,.0f} rec/s)")


if __name__ == "__main__":
    main()
