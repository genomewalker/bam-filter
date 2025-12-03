#!/usr/bin/env python3
"""
Simplest parallel SAM/BAM to Parquet converter.

Strategy: Use GNU parallel with samtools to split by reference,
then convert each chunk in parallel.
"""

import sys
import time
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed


def get_references(sam_file: str):
    """Get list of references from header."""
    result = subprocess.run(
        ['samtools', 'view', '-H', sam_file],
        capture_output=True,
        text=True,
        check=True
    )

    refs = []
    for line in result.stdout.splitlines():
        if line.startswith('@SQ'):
            for field in line.split('\t'):
                if field.startswith('SN:'):
                    refs.append(field[3:])
                    break

    return refs


def process_reference_chunk(args):
    """
    Process a chunk of references.

    Returns stats dict.
    """
    chunk_id, refs, input_file, output_dir, num_partitions, batch_size, compression_level, num_threads = args

    from bam_filter.parquet_converter_pure_cpp import convert_sam_bam_to_parquet_pure_cpp

    start_time = time.time()

    # Create temp filtered file using samtools
    temp_dir = Path(output_dir) / f"chunk_{chunk_id:03d}"
    temp_dir.mkdir(parents=True, exist_ok=True)
    temp_sam_gz = temp_dir / f"chunk_{chunk_id:03d}.sam.gz"

    print(f"[Chunk {chunk_id}] Filtering {len(refs):,} references...")

    # Get header first
    header_result = subprocess.run(
        ['samtools', 'view', '-H', input_file],
        capture_output=True,
        check=True
    )

    # Create set for fast lookup
    refs_set = set(refs)

    import gzip

    # Filter records: read full file, only keep assigned refs
    with gzip.open(temp_sam_gz, 'wt') as out:
        # Write header
        out.write(header_result.stdout.decode())

        # Stream and filter records
        cmd = ['samtools', 'view', input_file]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

        count = 0
        for line in proc.stdout:
            fields = line.split('\t', 4)  # Only split first 4 fields (faster)
            if len(fields) >= 3 and fields[2] in refs_set:
                out.write(line)
                count += 1
                if count % 500000 == 0:
                    print(f"[Chunk {chunk_id}] Filtered {count:,} records...")

        proc.wait()

    print(f"[Chunk {chunk_id}] Filtered {count:,} total records")

    filter_time = time.time() - start_time
    print(f"[Chunk {chunk_id}] Filtered in {filter_time:.1f}s, converting to Parquet...")

    # Convert to Parquet
    convert_start = time.time()
    stats = convert_sam_bam_to_parquet_pure_cpp(
        str(temp_sam_gz),
        str(temp_dir),
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
    temp_sam_gz.unlink()

    print(f"[Chunk {chunk_id}] DONE: {stats['total_records']:,} records in {total_time:.1f}s")

    return {
        'chunk_id': chunk_id,
        'num_refs': len(refs),
        'total_records': stats['total_records'],
        'total_time': total_time,
        'filter_time': filter_time,
        'convert_time': convert_time,
    }


def parallel_convert(
    input_file: str,
    output_dir: str,
    num_processes: int = 32,
    refs_per_chunk: int = 100000,  # 100K refs per chunk
    num_partitions: int = 16,
    batch_size: int = 100000,
    compression_level: int = 3,
    num_threads_per_worker: int = 2,
):
    """
    Convert SAM/BAM to Parquet in parallel by splitting references.

    Args:
        input_file: Input SAM/BAM/SAM.gz file
        output_dir: Output directory
        num_processes: Number of parallel processes
        refs_per_chunk: Number of references per chunk
        num_partitions: Partitions per chunk
        batch_size: Records per batch
        compression_level: ZSTD compression level
        num_threads_per_worker: Threads per worker

    Returns:
        Statistics dict
    """
    overall_start = time.time()

    print(f"Parallel SAM/BAM → Parquet Converter")
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print(f"Max processes: {num_processes}")
    print(f"Refs per chunk: {refs_per_chunk:,}")
    print()

    # Get references
    print("Reading references from header...")
    refs = get_references(input_file)
    print(f"Found {len(refs):,} references")

    # Split into chunks
    chunks = []
    for i in range(0, len(refs), refs_per_chunk):
        chunk_refs = refs[i:i + refs_per_chunk]
        chunks.append((len(chunks), chunk_refs, input_file, output_dir,
                      num_partitions, batch_size, compression_level, num_threads_per_worker))

    print(f"Split into {len(chunks)} chunks")
    print()

    # Process in parallel
    print(f"Processing {len(chunks)} chunks with up to {num_processes} parallel workers...")
    results = []

    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = {executor.submit(process_reference_chunk, chunk): chunk[0] for chunk in chunks}

        for future in as_completed(futures):
            chunk_id = futures[future]
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                print(f"[Chunk {chunk_id}] FAILED: {e}")

    total_time = time.time() - overall_start
    total_records = sum(r['total_records'] for r in results)

    print()
    print("="*80)
    print("PARALLEL CONVERSION COMPLETE")
    print("="*80)
    print(f"Total records: {total_records:,}")
    print(f"Total time: {total_time:.1f}s ({total_time/60:.1f} minutes)")
    print(f"Throughput: {total_records/total_time:.0f} records/sec")
    print(f"Throughput: {total_records/total_time/1e6:.2f} M records/sec")
    print()

    return {
        'total_records': total_records,
        'total_time_seconds': total_time,
        'throughput_records_per_sec': total_records / total_time,
        'num_chunks': len(chunks),
        'chunk_results': results,
    }


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python parallel_parquet_simple.py <input.sam.gz> <output_dir> [num_processes]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    num_processes = int(sys.argv[3]) if len(sys.argv) > 3 else 32

    stats = parallel_convert(
        input_file,
        output_dir,
        num_processes=num_processes,
        refs_per_chunk=200000,  # 200K refs per chunk
        num_partitions=16,
        batch_size=100000,
        compression_level=3,
        num_threads_per_worker=2,
    )
