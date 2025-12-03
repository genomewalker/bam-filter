#!/usr/bin/env python3
"""
Parallel SAM/BAM to Parquet converter with proper bgzip block handling.

Strategy: Split bgzip file at block boundaries (64KB blocks).
Each worker processes complete, valid bgzip blocks.
"""

import sys
import os
import struct
import time
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed


def find_bgzip_blocks(filepath, sample_rate=1000):
    """
    Scan bgzip file to find block offsets.

    BGZF format: each block has a gzip header with extra fields.
    We scan for block boundaries by reading the compressed size from each header.

    Returns list of (block_offset, block_size) tuples.
    """
    blocks = []
    offset = 0

    with open(filepath, 'rb') as f:
        block_count = 0
        while True:
            # Read gzip header (10 bytes minimum)
            header = f.read(10)
            if len(header) < 10:
                break

            # Check gzip magic number
            if header[0:2] != b'\x1f\x8b':
                print(f"Warning: Invalid gzip header at offset {offset}")
                break

            # Read extra field length (bytes 10-11)
            xlen = struct.unpack('<H', header[10-2:10])[0] if len(header) >= 12 else 0

            # Read extra fields to find BSIZE
            f.seek(offset + 10)
            extra = f.read(2)  # XLEN
            if len(extra) < 2:
                break
            xlen = struct.unpack('<H', extra)[0]

            # Read extra subfields to find BC (bgzip) subfield
            extra_data = f.read(xlen)

            # Parse extra subfields to find BSIZE
            # BC subfield: SI1='B', SI2='C', SLEN=2, BSIZE=uint16
            bsize = None
            i = 0
            while i < len(extra_data) - 4:
                si1 = extra_data[i]
                si2 = extra_data[i+1]
                slen = struct.unpack('<H', extra_data[i+2:i+4])[0]

                if si1 == ord('B') and si2 == ord('C'):
                    # Found BC subfield - contains BSIZE
                    if i + 4 + slen <= len(extra_data):
                        bsize = struct.unpack('<H', extra_data[i+4:i+6])[0]
                    break

                i += 4 + slen

            if bsize is None:
                # Not a bgzip block, might be regular gzip
                print(f"Warning: No BGZF BSIZE found at offset {offset}")
                break

            # BSIZE is the total block size minus 1
            block_size = bsize + 1

            # Store block info (sample to avoid memory issues)
            if block_count % sample_rate == 0:
                blocks.append((offset, block_size))

            block_count += 1

            # Skip to next block
            offset += block_size
            f.seek(offset)

            if block_count % 10000 == 0:
                print(f"  Scanned {block_count:,} blocks ({offset/1e6:.1f} MB)...")

    # Always include the last block
    blocks.append((offset, 0))  # End marker

    print(f"Found {block_count:,} total blocks, sampled {len(blocks):,} boundaries")
    return blocks, block_count


def process_bgzip_chunk(args):
    """Process a chunk of bgzip blocks."""
    chunk_id, input_file, start_offset, end_offset, output_dir, num_partitions, batch_size, compression_level, num_threads = args

    from bam_filter.parquet_converter_pure_cpp import convert_sam_bam_to_parquet_pure_cpp

    start_time = time.time()

    chunk_dir = Path(output_dir) / f"chunk_{chunk_id:03d}"
    chunk_dir.mkdir(parents=True, exist_ok=True)
    chunk_file = chunk_dir / f"chunk_{chunk_id:03d}.sam.gz"

    chunk_size = end_offset - start_offset
    chunk_size_mb = chunk_size / 1e6

    print(f"[Chunk {chunk_id}] Extracting bgzip blocks: offset {start_offset:,} to {end_offset:,} ({chunk_size_mb:.1f} MB)")

    # Extract byte range using dd at block boundaries
    with open(chunk_file, 'wb') as out:
        cmd = [
            'dd',
            f'if={input_file}',
            f'bs=1',
            f'skip={start_offset}',
            f'count={chunk_size}',
            'status=none'
        ]
        subprocess.run(cmd, stdout=out, check=True)

    extract_time = time.time() - start_time
    actual_size_mb = chunk_file.stat().st_size / 1e6
    print(f"[Chunk {chunk_id}] Extracted {actual_size_mb:.1f} MB in {extract_time:.1f}s")

    # Verify it's valid bgzip
    try:
        test_result = subprocess.run(
            ['bgzip', '-t', str(chunk_file)],
            capture_output=True,
            timeout=10
        )
        if test_result.returncode != 0:
            print(f"[Chunk {chunk_id}] WARNING: bgzip test failed - may have corrupt blocks")
    except Exception as e:
        print(f"[Chunk {chunk_id}] Could not verify bgzip: {e}")

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

        print(f"[Chunk {chunk_id}] DONE: {stats['total_records']:,} records in {total_time:.1f}s total ({convert_time:.1f}s convert)")

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
        import traceback
        traceback.print_exc()
        return {
            'chunk_id': chunk_id,
            'total_records': 0,
            'total_time': time.time() - start_time,
            'extract_time': extract_time,
            'convert_time': 0,
            'bytes_processed': chunk_size,
            'error': str(e),
        }


def parallel_convert_bgzip(
    input_file: str,
    output_dir: str,
    num_processes: int = 32,
    num_partitions: int = 16,
    batch_size: int = 100000,
    compression_level: int = 3,
    num_threads_per_worker: int = 2,
    block_sample_rate: int = 1000,
):
    """
    Convert bgzip SAM to Parquet in parallel, respecting bgzip block boundaries.

    Args:
        input_file: Input .sam.gz file (bgzip compressed)
        output_dir: Output directory
        num_processes: Number of parallel processes
        num_partitions: Partitions per chunk
        batch_size: Records per batch
        compression_level: ZSTD compression level
        num_threads_per_worker: Threads per worker
        block_sample_rate: Sample every Nth block (1000 = sample ~every 64MB)
    """
    overall_start = time.time()

    print(f"Parallel bgzip SAM → Parquet Converter (BGZF Block-Aware)")
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print(f"Processes: {num_processes}")
    print()

    # Scan bgzip blocks
    print("Scanning bgzip block boundaries...")
    scan_start = time.time()
    blocks, total_blocks = find_bgzip_blocks(input_file, sample_rate=block_sample_rate)
    scan_time = time.time() - scan_start

    file_size = os.path.getsize(input_file)
    print(f"File size: {file_size/1e9:.2f} GB ({file_size:,} bytes)")
    print(f"Total blocks: {total_blocks:,}")
    print(f"Sampled boundaries: {len(blocks):,}")
    print(f"Scan time: {scan_time:.1f}s")
    print()

    # Split sampled blocks among workers
    blocks_per_worker = len(blocks) // num_processes

    chunks = []
    for i in range(num_processes):
        start_idx = i * blocks_per_worker
        end_idx = start_idx + blocks_per_worker if i < num_processes - 1 else len(blocks) - 1

        start_offset = blocks[start_idx][0]
        end_offset = blocks[end_idx][0]

        chunks.append((
            i, input_file, start_offset, end_offset, output_dir,
            num_partitions, batch_size, compression_level, num_threads_per_worker
        ))

    print(f"Created {len(chunks)} chunks at bgzip block boundaries")
    for i, (chunk_id, _, start, end, *_) in enumerate(chunks[:5]):
        print(f"  Chunk {chunk_id}: bytes {start:,} to {end:,} ({(end-start)/1e6:.1f} MB)")
    if len(chunks) > 5:
        print(f"  ... and {len(chunks) - 5} more chunks")
    print()

    # Process in parallel
    print(f"Processing {len(chunks)} chunks with {num_processes} parallel workers...")
    results = []

    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = {executor.submit(process_bgzip_chunk, chunk): chunk[0] for chunk in chunks}

        for future in as_completed(futures):
            chunk_id = futures[future]
            try:
                result = future.result()
                results.append(result)
                print(f"✓ Chunk {chunk_id} completed: {result['total_records']:,} records")
            except Exception as e:
                print(f"✗ Chunk {chunk_id} EXCEPTION: {e}")

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
    print(f"Data rate: {file_size/1e9/total_time:.2f} GB/sec")
    print()

    return {
        'total_records': total_records,
        'total_time_seconds': total_time,
        'throughput_records_per_sec': total_records / total_time,
        'num_chunks': len(chunks),
        'chunk_results': results,
        'file_size_bytes': file_size,
        'total_blocks': total_blocks,
    }


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python parallel_parquet_bgzip.py <input.sam.gz> <output_dir> [num_processes]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    num_processes = int(sys.argv[3]) if len(sys.argv) > 3 else 32

    stats = parallel_convert_bgzip(
        input_file,
        output_dir,
        num_processes=num_processes,
        num_partitions=16,
        batch_size=100000,
        compression_level=3,
        num_threads_per_worker=2,
        block_sample_rate=1000,  # Sample every 1000 blocks (~64MB)
    )
