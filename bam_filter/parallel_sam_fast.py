#!/usr/bin/env python3
"""
Parallel SAM.gz to Parquet converter using the fast Cython parser.

Uses BGZF block-based splitting for true parallel processing.
Each worker uses the fast Cython parser that bypasses HTSlib.
"""

import os
import struct
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import time


def find_bgzf_blocks(filepath: str, num_chunks: int = 8) -> list:
    """
    Find BGZF block boundaries for parallel processing.

    Returns list of (start_offset, end_offset) tuples.
    """
    file_size = os.path.getsize(filepath)
    chunk_size = file_size // num_chunks

    boundaries = [0]

    with open(filepath, 'rb') as f:
        for i in range(1, num_chunks):
            target = i * chunk_size
            f.seek(target)

            # Search for BGZF block header (magic: 1f 8b 08 04)
            search_buffer = f.read(65536 * 2)

            pos = 0
            while pos < len(search_buffer) - 18:
                if search_buffer[pos:pos+4] == b'\x1f\x8b\x08\x04':
                    # Verify BC subfield
                    if pos + 12 < len(search_buffer):
                        xlen = struct.unpack('<H', search_buffer[pos+10:pos+12])[0]
                        if xlen >= 6 and pos + 12 + xlen <= len(search_buffer):
                            extra = search_buffer[pos+12:pos+12+xlen]
                            if extra[:2] == b'BC' and len(extra) >= 6:
                                boundaries.append(target + pos)
                                break
                pos += 1
            else:
                boundaries.append(target)

    boundaries.append(file_size)

    return [(boundaries[i], boundaries[i+1]) for i in range(len(boundaries)-1)]


def process_sam_range(args):
    """Process a BGZF byte range using the fast Cython parser."""
    (
        worker_id, input_file, start_offset, end_offset, output_dir,
        batch_size, compression_level, store_sequences
    ) = args

    from bam_filter.sam_parser_fast import convert_sam_gz_range_to_parquet_fast

    start_time = time.time()
    worker_dir = Path(output_dir) / f"worker_{worker_id:03d}"
    worker_dir.mkdir(parents=True, exist_ok=True)

    print(f"[Worker {worker_id}] Processing bytes {start_offset:,} to {end_offset:,}")

    try:
        stats = convert_sam_gz_range_to_parquet_fast(
            input_file,
            str(worker_dir),
            start_offset,
            end_offset,
            batch_size=batch_size,
            compression_level=compression_level,
            num_threads=1,  # Each worker uses 1 thread
            store_sequences=store_sequences,
        )

        elapsed = time.time() - start_time
        records = stats.get('total_records', 0)
        rate = records / elapsed if elapsed > 0 else 0

        print(f"[Worker {worker_id}] Done: {records:,} records in {elapsed:.1f}s ({rate:,.0f} rec/s)")

        return {
            'worker_id': worker_id,
            'total_records': records,
            'skipped_lines': stats.get('skipped_lines', 0),
            'elapsed': elapsed,
            'output_dir': str(worker_dir),
        }
    except Exception as e:
        import traceback
        print(f"[Worker {worker_id}] Error: {e}")
        traceback.print_exc()
        return {
            'worker_id': worker_id,
            'total_records': 0,
            'elapsed': time.time() - start_time,
            'error': str(e),
        }


def parallel_convert_sam_fast(
    input_file: str,
    output_dir: str,
    num_workers: int = 8,
    batch_size: int = 100000,
    compression_level: int = 3,
    store_sequences: bool = False,
) -> dict:
    """
    Convert SAM.gz to Parquet using parallel fast Cython parser.

    Each worker:
    - Seeks to its BGZF block range
    - Uses the fast Cython parser (bypasses HTSlib)
    - Writes to its own Parquet file

    Args:
        input_file: Path to SAM.gz file
        output_dir: Output directory
        num_workers: Number of parallel workers
        batch_size: Records per batch
        compression_level: ZSTD compression level
        store_sequences: Store sequence/quality strings

    Returns:
        Statistics dict
    """
    start_time = time.time()

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Find BGZF block boundaries
    print(f"Scanning BGZF blocks for {num_workers} workers...")
    ranges = find_bgzf_blocks(input_file, num_workers)
    print(f"  Found {len(ranges)} ranges")

    # Prepare worker arguments
    worker_args = [
        (
            i, input_file, start, end, output_dir,
            batch_size, compression_level, store_sequences
        )
        for i, (start, end) in enumerate(ranges)
    ]

    # Process in parallel
    print(f"\nStarting {num_workers} workers with fast Cython parser...")
    results = []

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(process_sam_range, args): args[0] for args in worker_args}

        for future in as_completed(futures):
            result = future.result()
            results.append(result)

    # Aggregate results
    total_records = sum(r.get('total_records', 0) for r in results)
    total_elapsed = time.time() - start_time

    print(f"\n{'='*60}")
    print(f"TOTAL: {total_records:,} records in {total_elapsed:.1f}s")
    print(f"Throughput: {total_records/total_elapsed:,.0f} rec/s ({total_records/total_elapsed/1e6:.3f} M/s)")
    print(f"{'='*60}")

    return {
        'total_records': total_records,
        'total_elapsed': total_elapsed,
        'throughput': total_records / total_elapsed if total_elapsed > 0 else 0,
        'worker_results': results,
    }


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 3:
        print("Usage: python parallel_sam_fast.py <input.sam.gz> <output_dir> [num_workers]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    num_workers = int(sys.argv[3]) if len(sys.argv) > 3 else 8

    stats = parallel_convert_sam_fast(
        input_file,
        output_dir,
        num_workers=num_workers,
        store_sequences=False,
    )
