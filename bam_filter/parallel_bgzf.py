#!/usr/bin/env python3
"""
Smart parallel BAM/SAM.gz to Parquet converter using BGZF block-based splitting.

Key Design Decisions:
1. SINGLE output partition (no redundant by_ref + by_read)
2. BGZF block-based parallel reading (seek directly, no record skipping)
3. Workers write to separate files, then merge
4. Minimal schema for fast processing

Architecture:
    Main Process
        ├── Scan BGZF blocks (find block boundaries)
        ├── Assign byte ranges to workers
        └── Workers (ProcessPoolExecutor)
            ├── Seek to block boundary
            ├── Read records until end boundary
            ├── Handle record overlap at boundaries
            └── Write to worker-specific Parquet file
"""

import os
import struct
import gzip
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import time


def find_bgzf_blocks(filepath: str, num_chunks: int = 8) -> list:
    """
    Find BGZF block boundaries for parallel processing.

    BGZF format: Each block starts with magic bytes 1f 8b 08 04
    followed by extra field with BC subfield containing block size.

    Returns list of (start_offset, end_offset) tuples for each chunk.
    """
    file_size = os.path.getsize(filepath)
    chunk_size = file_size // num_chunks

    boundaries = [0]

    with open(filepath, 'rb') as f:
        for i in range(1, num_chunks):
            # Seek to approximate chunk boundary
            target = i * chunk_size
            f.seek(target)

            # Search for next BGZF block header (magic: 1f 8b 08 04)
            # Read in chunks to find the pattern
            search_buffer = f.read(65536 * 2)  # Search up to 128KB ahead

            # Look for BGZF magic bytes
            pos = 0
            while pos < len(search_buffer) - 18:
                if (search_buffer[pos:pos+4] == b'\x1f\x8b\x08\x04'):
                    # Verify it's a valid BGZF block by checking extra field
                    # Skip to extra field length at offset 10
                    if pos + 12 < len(search_buffer):
                        xlen = struct.unpack('<H', search_buffer[pos+10:pos+12])[0]
                        if xlen >= 6 and pos + 12 + xlen <= len(search_buffer):
                            # Check for BC subfield
                            extra = search_buffer[pos+12:pos+12+xlen]
                            if extra[:2] == b'BC' and len(extra) >= 6:
                                # Valid BGZF block found
                                boundaries.append(target + pos)
                                break
                pos += 1
            else:
                # No block found, use approximate boundary
                boundaries.append(target)

    boundaries.append(file_size)

    # Create (start, end) ranges
    ranges = [(boundaries[i], boundaries[i+1]) for i in range(len(boundaries)-1)]
    return ranges


def process_bgzf_range(args):
    """Process a BGZF byte range in parallel - calls optimized Cython code."""
    (
        worker_id, input_file, start_offset, end_offset, output_dir,
        batch_size, compression_level, num_threads, calculate_pmd,
        library_type, store_sequences
    ) = args

    from bam_filter.parquet_converter_pure_cpp import convert_bgzf_range_to_parquet

    start_time = time.time()
    worker_dir = Path(output_dir) / f"worker_{worker_id:03d}"
    worker_dir.mkdir(parents=True, exist_ok=True)

    print(f"[Worker {worker_id}] Processing bytes {start_offset:,} to {end_offset:,}")

    try:
        stats = convert_bgzf_range_to_parquet(
            input_file,
            str(worker_dir),
            start_offset,
            end_offset,
            batch_size=batch_size,
            compression_level=compression_level,
            num_threads=num_threads,
            calculate_pmd=calculate_pmd,
            library_type=library_type,
            store_sequences=store_sequences,
        )

        elapsed = time.time() - start_time
        records = stats.get('total_records', 0)
        rate = records / elapsed if elapsed > 0 else 0

        print(f"[Worker {worker_id}] Done: {records:,} records in {elapsed:.1f}s ({rate:,.0f} rec/s)")

        return {
            'worker_id': worker_id,
            'total_records': records,
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


def parallel_convert_bgzf(
    input_file: str,
    output_dir: str,
    num_workers: int = 8,
    batch_size: int = 100000,
    compression_level: int = 3,
    num_threads_per_worker: int = 1,
    calculate_pmd: bool = False,
    library_type: str = "ds",
    store_sequences: bool = False,
) -> dict:
    """
    Convert BAM/SAM.gz to Parquet using parallel BGZF block processing.

    This is the FAST approach:
    - Each worker seeks directly to its byte range (no record skipping!)
    - Workers process in parallel without contention
    - Results are written to separate files, then can be merged

    Args:
        input_file: Path to BAM or SAM.gz file
        output_dir: Output directory for Parquet files
        num_workers: Number of parallel workers
        batch_size: Records per batch
        compression_level: ZSTD compression level (1-22)
        num_threads_per_worker: HTSlib threads per worker
        calculate_pmd: Calculate PMD scores
        library_type: Library type for PMD ("ds" or "ss")
        store_sequences: Store sequence/quality strings

    Returns:
        Statistics dict with total records and timing
    """
    from bam_filter.parquet_converter_pure_cpp import save_bam_header

    start_time = time.time()

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Save header
    print(f"Saving header from {input_file}...")
    save_bam_header(input_file, output_dir)

    # Find BGZF block boundaries
    print(f"Scanning BGZF blocks for {num_workers} workers...")
    ranges = find_bgzf_blocks(input_file, num_workers)
    print(f"  Found {len(ranges)} ranges")

    # Prepare worker arguments
    worker_args = [
        (
            i, input_file, start, end, output_dir,
            batch_size, compression_level, num_threads_per_worker,
            calculate_pmd, library_type, store_sequences
        )
        for i, (start, end) in enumerate(ranges)
    ]

    # Process in parallel
    print(f"\nStarting {num_workers} workers...")
    results = []

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(process_bgzf_range, args): args[0] for args in worker_args}

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
        print("Usage: python parallel_bgzf.py <input.bam|sam.gz> <output_dir> [num_workers]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    num_workers = int(sys.argv[3]) if len(sys.argv) > 3 else 8

    stats = parallel_convert_bgzf(
        input_file,
        output_dir,
        num_workers=num_workers,
        calculate_pmd=False,
        store_sequences=False,
    )
