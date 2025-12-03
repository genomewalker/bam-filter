#!/usr/bin/env python3
"""
Parallel SAM/BAM to Parquet converter using byte-range splitting.

Strategy: Split the bgzip file by byte ranges, each worker processes its chunk.
This avoids reading the file 32x times.
"""

import sys
import os
import time
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed


def get_file_size(filepath):
    """Get file size in bytes."""
    return os.path.getsize(filepath)


def process_byte_chunk(args):
    """
    Process a byte range chunk of the SAM/BAM file.

    Returns stats dict.
    """
    chunk_id, input_file, start_byte, end_byte, output_dir, num_partitions, batch_size, compression_level, num_threads = args

    from bam_filter.parquet_converter_pure_cpp import convert_sam_bam_to_parquet_pure_cpp

    start_time = time.time()

    # Create chunk-specific directory
    chunk_dir = Path(output_dir) / f"chunk_{chunk_id:03d}"
    chunk_dir.mkdir(parents=True, exist_ok=True)

    chunk_file = chunk_dir / f"chunk_{chunk_id:03d}.sam.gz"

    print(f"[Chunk {chunk_id}] Extracting bytes {start_byte:,} to {end_byte:,} ({(end_byte-start_byte)/1e6:.1f} MB)...")

    # Extract byte range using dd
    # dd if=input.sam.gz of=chunk.sam.gz bs=1M skip=X count=Y
    chunk_size = end_byte - start_byte

    # Use dd to extract byte range
    with open(chunk_file, 'wb') as out:
        cmd = [
            'dd',
            f'if={input_file}',
            'bs=1M',
            f'skip={start_byte // (1024*1024)}',  # Skip in MB
            f'count={(chunk_size // (1024*1024)) + 1}',  # Count in MB (add 1 to ensure we get everything)
            'status=none'
        ]
        subprocess.run(cmd, stdout=out, check=True)

    extract_time = time.time() - start_time
    file_size_mb = chunk_file.stat().st_size / 1e6
    print(f"[Chunk {chunk_id}] Extracted {file_size_mb:.1f} MB in {extract_time:.1f}s, converting...")

    # Convert to Parquet
    convert_start = time.time()

    try:
        stats = convert_sam_bam_to_parquet_pure_cpp(
            str(chunk_file),
            str(chunk_dir),
            num_partitions=num_partitions,
            batch_size=batch_size,
            write_by_reference=True,
            write_by_read=True,
            compression_level=compression_level,
            num_threads=num_threads,
        )

        convert_time = time.time() - convert_start
        total_time = time.time() - start_time

        # Clean up temp file
        chunk_file.unlink()

        print(f"[Chunk {chunk_id}] DONE: {stats['total_records']:,} records in {total_time:.1f}s")

        return {
            'chunk_id': chunk_id,
            'total_records': stats['total_records'],
            'total_time': total_time,
            'extract_time': extract_time,
            'convert_time': convert_time,
            'bytes_processed': chunk_size,
        }
    except Exception as e:
        print(f"[Chunk {chunk_id}] FAILED: {e}")
        # Return zero records on failure
        return {
            'chunk_id': chunk_id,
            'total_records': 0,
            'total_time': time.time() - start_time,
            'extract_time': extract_time,
            'convert_time': 0,
            'bytes_processed': chunk_size,
            'error': str(e),
        }


def parallel_convert_byterange(
    input_file: str,
    output_dir: str,
    num_processes: int = 32,
    num_partitions: int = 16,
    batch_size: int = 100000,
    compression_level: int = 3,
    num_threads_per_worker: int = 2,
):
    """
    Convert SAM/BAM to Parquet in parallel by splitting file into byte ranges.

    Args:
        input_file: Input SAM/BAM/SAM.gz file
        output_dir: Output directory
        num_processes: Number of parallel processes
        num_partitions: Partitions per chunk
        batch_size: Records per batch
        compression_level: ZSTD compression level
        num_threads_per_worker: Threads per worker

    Returns:
        Statistics dict
    """
    overall_start = time.time()

    print(f"Parallel SAM/BAM → Parquet Converter (Byte-Range Splitting)")
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print(f"Processes: {num_processes}")
    print()

    # Get file size
    file_size = get_file_size(input_file)
    file_size_mb = file_size / 1e6
    file_size_gb = file_size / 1e9

    print(f"File size: {file_size_gb:.2f} GB ({file_size:,} bytes)")

    # Calculate chunk size
    chunk_size = file_size // num_processes
    chunk_size_mb = chunk_size / 1e6

    print(f"Chunk size: {chunk_size_mb:.1f} MB ({chunk_size:,} bytes)")
    print()

    # Create byte range chunks
    chunks = []
    for i in range(num_processes):
        start_byte = i * chunk_size
        end_byte = start_byte + chunk_size if i < num_processes - 1 else file_size

        chunks.append((
            i, input_file, start_byte, end_byte, output_dir,
            num_partitions, batch_size, compression_level, num_threads_per_worker
        ))

    print(f"Created {len(chunks)} byte-range chunks")
    print()

    # Process in parallel
    print(f"Processing {len(chunks)} chunks with {num_processes} parallel workers...")
    results = []

    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = {executor.submit(process_byte_chunk, chunk): chunk[0] for chunk in chunks}

        for future in as_completed(futures):
            chunk_id = futures[future]
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                print(f"[Chunk {chunk_id}] EXCEPTION: {e}")

    total_time = time.time() - overall_start
    total_records = sum(r['total_records'] for r in results)

    print()
    print("="*80)
    print("PARALLEL CONVERSION COMPLETE")
    print("="*80)
    print(f"Total records: {total_records:,}")
    print(f"Total time: {total_time:.1f}s ({total_time/60:.1f} minutes)")
    print(f"Throughput: {total_records/total_time:.0f} records/sec")
    print(f"Throughput: {total_records/total_time/1e6:.3f} M records/sec")
    print(f"Data rate: {file_size_gb/total_time:.1f} GB/sec")
    print()

    return {
        'total_records': total_records,
        'total_time_seconds': total_time,
        'throughput_records_per_sec': total_records / total_time,
        'num_chunks': len(chunks),
        'chunk_results': results,
        'file_size_bytes': file_size,
    }


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python parallel_parquet_byterange.py <input.sam.gz> <output_dir> [num_processes]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    num_processes = int(sys.argv[3]) if len(sys.argv) > 3 else 32

    stats = parallel_convert_byterange(
        input_file,
        output_dir,
        num_processes=num_processes,
        num_partitions=16,
        batch_size=100000,
        compression_level=3,
        num_threads_per_worker=2,
    )
