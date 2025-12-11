#!/usr/bin/env python3
"""
Test the fast Cython SAM.gz parser that bypasses HTSlib text parsing.

Usage:
    python test_sam_parser_fast.py <input.sam.gz> <output_dir>
"""

import sys
import time
from pathlib import Path


def main():
    if len(sys.argv) < 3:
        print("Usage: python test_sam_parser_fast.py <input.sam.gz> <output_dir>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]

    from bam_filter.sam_parser_fast import convert_sam_gz_to_parquet_fast

    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print(f"Options: Fast Cython parser, NO SEQUENCES, 8 threads for BGZF")
    print()

    # Clean output
    output_path = Path(output_dir)
    if output_path.exists():
        import shutil
        shutil.rmtree(output_path)

    start_time = time.time()

    stats = convert_sam_gz_to_parquet_fast(
        input_file,
        output_dir,
        batch_size=100000,
        compression_level=3,
        num_threads=8,
        store_sequences=False,
    )

    total_elapsed = time.time() - start_time
    total_records = stats['total_records']

    print(f"\n{'='*60}")
    print(f"FAST SAM PARSER RESULTS")
    print(f"{'='*60}")
    print(f"Total records: {total_records:,}")
    print(f"Total time: {total_elapsed:.1f}s")
    print(f"Throughput: {total_records/total_elapsed:,.0f} rec/s")
    print(f"           {total_records/total_elapsed/1e6:.3f} M/s")

    # Output size
    parquet_files = list(Path(output_dir).rglob("*.parquet"))
    total_size = sum(f.stat().st_size for f in parquet_files)
    print(f"Output: {len(parquet_files)} files, {total_size/1e6:.1f} MB")


if __name__ == "__main__":
    main()
