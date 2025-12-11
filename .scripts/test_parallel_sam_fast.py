#!/usr/bin/env python3
"""
Test parallel SAM.gz to Parquet using the fast Cython parser.

Usage:
    python test_parallel_sam_fast.py <input.sam.gz> <output_dir> [num_workers]
"""

import sys
import time
from pathlib import Path


def main():
    if len(sys.argv) < 3:
        print("Usage: python test_parallel_sam_fast.py <input.sam.gz> <output_dir> [num_workers]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    num_workers = int(sys.argv[3]) if len(sys.argv) > 3 else 8

    from bam_filter.parallel_sam_fast import parallel_convert_sam_fast

    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print(f"Workers: {num_workers}")
    print(f"Options: Fast Cython parser (bypasses HTSlib), NO SEQUENCES")
    print()

    # Clean output
    output_path = Path(output_dir)
    if output_path.exists():
        import shutil
        shutil.rmtree(output_path)

    start_time = time.time()

    stats = parallel_convert_sam_fast(
        input_file,
        output_dir,
        num_workers=num_workers,
        batch_size=100000,
        compression_level=3,
        store_sequences=False,
    )

    total_elapsed = time.time() - start_time
    total_records = stats['total_records']

    print(f"\n{'='*60}")
    print(f"PARALLEL FAST SAM PARSER RESULTS")
    print(f"{'='*60}")
    print(f"Total records: {total_records:,}")
    print(f"Total time: {total_elapsed:.1f}s")
    print(f"Throughput: {total_records/total_elapsed:,.0f} rec/s")
    print(f"           {total_records/total_elapsed/1e6:.3f} M/s")

    # Output size
    parquet_files = list(Path(output_dir).rglob("*.parquet"))
    total_size = sum(f.stat().st_size for f in parquet_files)
    print(f"Output: {len(parquet_files)} files, {total_size/1e6:.1f} MB")

    # Per-worker stats
    print(f"\nPer-worker breakdown:")
    for r in sorted(stats['worker_results'], key=lambda x: x['worker_id']):
        worker_id = r['worker_id']
        records = r['total_records']
        elapsed = r['elapsed']
        rate = records / elapsed if elapsed > 0 else 0
        skipped = r.get('skipped_lines', 0)
        if 'error' in r:
            print(f"  Worker {worker_id}: ERROR - {r['error']}")
        else:
            print(f"  Worker {worker_id}: {records:>12,} records in {elapsed:>6.1f}s ({rate:>9,.0f} rec/s) [skipped: {skipped}]")


if __name__ == "__main__":
    main()
