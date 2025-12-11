#!/usr/bin/env python3
"""Benchmark BAM to Parquet - NO SEQUENCES mode for maximum speed."""

import sys
import time
from pathlib import Path

def main():
    if len(sys.argv) < 3:
        print("Usage: python test_bam_to_parquet_noseq.py <input.bam> <output_dir>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]

    from bam_filter.parquet_converter_pure_cpp import (
        convert_sam_bam_to_parquet_pure_cpp,
        save_bam_header,
    )

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print(f"Options: NO PMD, NO SEQUENCES, by_reference only, 8 threads")
    print()

    save_bam_header(input_file, output_dir)

    start_time = time.time()

    stats = convert_sam_bam_to_parquet_pure_cpp(
        input_file,
        output_dir,
        num_partitions=16,
        batch_size=100000,
        write_by_reference=True,
        write_by_read=False,
        compression_level=3,
        num_threads=8,
        skip_records=0,
        max_records=0,
        calculate_pmd=False,
        library_type="ds",
        store_sequences=False,  # Skip sequence/quality storage
    )

    elapsed = time.time() - start_time
    records = stats["total_records"]
    rate = records / elapsed if elapsed > 0 else 0

    print(f"\nRecords: {records:,}")
    print(f"Time: {elapsed:.1f}s")
    print(f"Throughput: {rate:,.0f} rec/s ({rate/1e6:.3f} M/s)")

    # Output size
    output_path = Path(output_dir)
    parquet_files = list(output_path.rglob("*.parquet"))
    total_size = sum(f.stat().st_size for f in parquet_files)
    print(f"Output: {len(parquet_files)} files, {total_size/1e6:.1f} MB")

if __name__ == "__main__":
    main()
