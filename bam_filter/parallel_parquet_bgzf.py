#!/usr/bin/env python3
"""
BGZF Block-based Parallel BAM/SAM to Parquet Converter.

This module uses BGZF block boundaries to correctly split compressed
BAM/SAM files across parallel workers. Each worker processes a range
of BGZF blocks, ensuring all records are captured.

Unlike the record-estimation approach, this method:
1. Scans BGZF block headers to find exact block boundaries
2. Assigns equal numbers of blocks to each worker
3. Workers seek to their starting block offset using HTSlib
4. Records spanning block boundaries are handled correctly
"""

import os
import sys
import time
import multiprocessing as mp
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any

from bam_filter.bgzf_scanner import (
    scan_bgzf_blocks,
    get_block_ranges,
    is_bgzf_file,
    print_bgzf_stats,
)


def process_bgzf_chunk(args: Tuple) -> Dict[str, Any]:
    """
    Process a chunk of BGZF blocks and convert to Parquet.

    Args:
        args: Tuple of (chunk_id, input_file, start_offset, end_offset,
              num_blocks, output_dir, num_partitions, batch_size,
              compression_level, num_threads, calculate_pmd, library_type)

    Returns:
        Dictionary with processing statistics
    """
    (chunk_id, input_file, start_offset, end_offset, num_blocks,
     output_dir, num_partitions, batch_size, compression_level,
     num_threads, calculate_pmd, library_type) = args

    # Import here to avoid issues with multiprocessing
    from bam_filter.parquet_converter_pure_cpp import convert_bgzf_range_to_parquet

    chunk_dir = Path(output_dir) / f"chunk_{chunk_id:03d}"
    chunk_dir.mkdir(parents=True, exist_ok=True)

    print(f"[Chunk {chunk_id}] Processing BGZF blocks at offset {start_offset:,} - {end_offset:,} ({num_blocks:,} blocks)")

    start_time = time.time()

    try:
        stats = convert_bgzf_range_to_parquet(
            input_file,
            str(chunk_dir),
            start_offset=start_offset,
            end_offset=end_offset,
            num_partitions=num_partitions,
            batch_size=batch_size,
            compression_level=compression_level,
            num_threads=num_threads,
            calculate_pmd=calculate_pmd,
            library_type=library_type,
        )

        elapsed = time.time() - start_time

        return {
            'chunk_id': chunk_id,
            'success': True,
            'records_processed': stats.get('records_processed', 0),
            'elapsed_seconds': elapsed,
            'start_offset': start_offset,
            'end_offset': end_offset,
            'num_blocks': num_blocks,
        }

    except Exception as e:
        elapsed = time.time() - start_time
        print(f"[Chunk {chunk_id}] ERROR: {e}")
        import traceback
        traceback.print_exc()

        return {
            'chunk_id': chunk_id,
            'success': False,
            'error': str(e),
            'elapsed_seconds': elapsed,
            'start_offset': start_offset,
            'end_offset': end_offset,
            'num_blocks': num_blocks,
        }


def parallel_convert_bgzf(
    input_file: str,
    output_dir: str,
    num_processes: int = 32,
    num_partitions: int = 16,
    batch_size: int = 100000,
    compression_level: int = 3,
    num_threads_per_worker: int = 1,
    calculate_pmd: bool = True,
    library_type: str = "ds",
) -> Dict[str, Any]:
    """
    Convert BGZF-compressed BAM/SAM to Parquet using block-based parallelization.

    This method:
    1. Scans the file to find all BGZF block boundaries
    2. Divides blocks equally among workers
    3. Each worker processes its assigned blocks
    4. Results are combined

    Args:
        input_file: Path to BGZF-compressed BAM/SAM file
        output_dir: Output directory for Parquet files
        num_processes: Number of parallel workers
        num_partitions: Partitions per chunk for Parquet
        batch_size: Records per batch
        compression_level: Parquet compression level
        num_threads_per_worker: Threads per worker (usually 1)
        calculate_pmd: Whether to calculate PMD scores
        library_type: "ds" (double-stranded) or "ss" (single-stranded)

    Returns:
        Dictionary with conversion statistics
    """
    print("="*80)
    print("BGZF Block-based Parallel BAM/SAM → Parquet Converter")
    print("="*80)
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print(f"Workers: {num_processes}")
    print()

    # Verify input is BGZF
    if not is_bgzf_file(input_file):
        raise ValueError(f"{input_file} is not a BGZF-compressed file")

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Save header first
    print("Saving BAM header for round-trip conversion...")
    from bam_filter.parquet_converter_pure_cpp import save_bam_header
    save_bam_header(input_file, output_dir)
    print()

    # Scan BGZF blocks
    print("Scanning BGZF blocks...")
    scan_start = time.time()
    block_offsets = scan_bgzf_blocks(input_file, progress_interval=1000000)
    scan_time = time.time() - scan_start
    print()

    # Print stats
    print_bgzf_stats(input_file, block_offsets)
    print()

    # Get file size
    file_size = Path(input_file).stat().st_size

    # Calculate block ranges for workers
    block_ranges = get_block_ranges(block_offsets, num_processes, file_size)

    print(f"Block distribution for {num_processes} workers:")
    for worker_id, start, end, blocks in block_ranges[:3]:
        print(f"  Worker {worker_id}: offset {start:,} - {end:,} ({blocks:,} blocks)")
    if len(block_ranges) > 3:
        print(f"  ... and {len(block_ranges) - 3} more workers")
    print()

    # Prepare worker arguments
    worker_args = []
    for worker_id, start_offset, end_offset, num_blocks in block_ranges:
        worker_args.append((
            worker_id,
            input_file,
            start_offset,
            end_offset,
            num_blocks,
            output_dir,
            num_partitions,
            batch_size,
            compression_level,
            num_threads_per_worker,
            calculate_pmd,
            library_type,
        ))

    # Process in parallel
    print(f"Processing with {len(worker_args)} parallel workers...")
    print("Each worker processes its assigned BGZF blocks")
    print()

    process_start = time.time()

    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(process_bgzf_chunk, worker_args)

    process_time = time.time() - process_start

    # Collect results
    total_records = 0
    successful_chunks = 0
    failed_chunks = []
    worker_times = []

    for result in results:
        if result['success']:
            total_records += result['records_processed']
            successful_chunks += 1
            worker_times.append(result['elapsed_seconds'])
            print(f"✓ Worker {result['chunk_id']} completed: {result['records_processed']:,} records")
        else:
            failed_chunks.append(result)
            print(f"✗ Worker {result['chunk_id']} FAILED: {result.get('error', 'unknown')}")

    print()
    print("="*80)
    print("BGZF PARALLEL CONVERSION COMPLETE")
    print("="*80)
    print(f"Total blocks scanned: {len(block_offsets):,}")
    print(f"Block scan time: {scan_time:.1f}s")
    print(f"Total records processed: {total_records:,}")
    print(f"Processing time: {process_time:.1f}s ({process_time/60:.2f} minutes)")
    print(f"Throughput: {total_records/process_time:,.0f} records/sec")
    print(f"Throughput: {total_records/process_time/1e6:.3f} M records/sec")
    print()

    if worker_times:
        avg_worker_time = sum(worker_times) / len(worker_times)
        print(f"Average worker time: {avg_worker_time:.1f}s")

    if failed_chunks:
        print(f"WARNING: {len(failed_chunks)} chunks failed!")
        for fc in failed_chunks:
            print(f"  Chunk {fc['chunk_id']}: {fc.get('error', 'unknown')}")

    return {
        'total_records': total_records,
        'total_blocks': len(block_offsets),
        'scan_time_seconds': scan_time,
        'process_time_seconds': process_time,
        'total_time_seconds': scan_time + process_time,
        'throughput_records_per_sec': total_records / process_time if process_time > 0 else 0,
        'num_workers': len(worker_args),
        'successful_chunks': successful_chunks,
        'failed_chunks': len(failed_chunks),
    }


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="BGZF Block-based Parallel Converter")
    parser.add_argument('input_file', help="Input BGZF file")
    parser.add_argument('output_dir', help="Output directory")
    parser.add_argument('-t', '--threads', type=int, default=32, help="Number of workers")
    parser.add_argument('--pmd', action='store_true', default=True, help="Calculate PMD")
    parser.add_argument('--library', choices=['ds', 'ss'], default='ds', help="Library type")

    args = parser.parse_args()

    stats = parallel_convert_bgzf(
        args.input_file,
        args.output_dir,
        num_processes=args.threads,
        calculate_pmd=args.pmd,
        library_type=args.library,
    )

    print()
    print(f"Conversion complete! {stats['total_records']:,} records")
