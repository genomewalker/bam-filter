#!/usr/bin/env python3
"""
BGZF (Blocked GNU Zip Format) block scanner for parallel BAM/SAM processing.

BGZF is the compression format used by BAM files and gzipped SAM files from
samtools. It consists of concatenated gzip blocks, each up to 64KB when
decompressed.

This module provides fast scanning of BGZF block boundaries, enabling
parallel processing by distributing blocks across workers.
"""

import struct
import time
from typing import List, Tuple, Optional
from pathlib import Path


# BGZF block header constants
GZIP_MAGIC = b'\x1f\x8b'
BGZF_EXTRA_FIELD_FLAG = 0x04
BGZF_SUBFIELD_ID = b'BC'


def scan_bgzf_blocks(
    file_path: str,
    progress_interval: int = 1000000,
    max_blocks: Optional[int] = None
) -> List[int]:
    """
    Scan a BGZF file and return all block start offsets.

    This is fast because it only reads block headers (18 bytes per block),
    not the actual compressed data.

    Args:
        file_path: Path to BGZF-compressed file
        progress_interval: Print progress every N blocks (0 to disable)
        max_blocks: Stop after scanning this many blocks (for testing)

    Returns:
        List of byte offsets where each BGZF block starts
    """
    block_offsets = []

    start_time = time.time()

    with open(file_path, 'rb') as f:
        # Get file size
        f.seek(0, 2)
        file_size = f.tell()
        f.seek(0)

        offset = 0
        count = 0

        while offset < file_size:
            if max_blocks and count >= max_blocks:
                break

            f.seek(offset)
            header = f.read(18)  # Minimum BGZF header size

            if len(header) < 18:
                break

            # Verify gzip magic bytes
            if header[0:2] != GZIP_MAGIC:
                raise ValueError(f"Invalid gzip magic at offset {offset}")

            # Check for extra field flag (bit 2 of FLG byte at position 3)
            if not (header[3] & BGZF_EXTRA_FIELD_FLAG):
                raise ValueError(f"No extra field flag at offset {offset}")

            # Verify BC subfield identifier
            if header[12:14] != BGZF_SUBFIELD_ID:
                raise ValueError(f"Not BGZF format (no BC subfield) at offset {offset}")

            # Get block size from BSIZE field (positions 16-17)
            # BSIZE is total block size - 1
            bsize = struct.unpack('<H', header[16:18])[0]
            block_size = bsize + 1

            block_offsets.append(offset)
            count += 1
            offset += block_size

            if progress_interval and count % progress_interval == 0:
                elapsed = time.time() - start_time
                print(f"  Scanned {count:,} blocks, {offset/1e9:.2f} GB, {elapsed:.1f}s")

    elapsed = time.time() - start_time
    if progress_interval:
        print(f"  Scan complete: {len(block_offsets):,} blocks in {elapsed:.1f}s")

    return block_offsets


def get_block_ranges(
    block_offsets: List[int],
    num_workers: int,
    file_size: int
) -> List[Tuple[int, int, int, int]]:
    """
    Divide BGZF blocks among workers.

    Args:
        block_offsets: List of block start offsets
        num_workers: Number of workers
        file_size: Total file size in bytes

    Returns:
        List of (worker_id, start_offset, end_offset, num_blocks) tuples
    """
    total_blocks = len(block_offsets)
    blocks_per_worker = total_blocks // num_workers
    remainder = total_blocks % num_workers

    ranges = []
    block_idx = 0

    for worker_id in range(num_workers):
        # Distribute remainder blocks to first workers
        worker_blocks = blocks_per_worker + (1 if worker_id < remainder else 0)

        if worker_blocks == 0:
            continue

        start_block = block_idx
        end_block = block_idx + worker_blocks

        start_offset = block_offsets[start_block]

        # End offset is start of next block or file size
        if end_block < len(block_offsets):
            end_offset = block_offsets[end_block]
        else:
            end_offset = file_size

        ranges.append((worker_id, start_offset, end_offset, worker_blocks))
        block_idx = end_block

    return ranges


def estimate_records_from_blocks(
    num_blocks: int,
    avg_records_per_block: float = 206.0
) -> int:
    """
    Estimate number of records based on block count.

    Each BGZF block decompresses to ~64KB. With an average SAM record
    size of ~300 bytes, that's roughly 206 records per block.

    This is a rough estimate - actual count will vary by file.

    Args:
        num_blocks: Number of BGZF blocks
        avg_records_per_block: Average records per block (default 206)

    Returns:
        Estimated record count
    """
    return int(num_blocks * avg_records_per_block)


def is_bgzf_file(file_path: str) -> bool:
    """
    Check if a file is BGZF-compressed.

    Args:
        file_path: Path to file

    Returns:
        True if file is BGZF format
    """
    try:
        with open(file_path, 'rb') as f:
            header = f.read(18)

            if len(header) < 18:
                return False

            # Check gzip magic
            if header[0:2] != GZIP_MAGIC:
                return False

            # Check extra field flag
            if not (header[3] & BGZF_EXTRA_FIELD_FLAG):
                return False

            # Check BC subfield
            if header[12:14] != BGZF_SUBFIELD_ID:
                return False

            return True

    except Exception:
        return False


def print_bgzf_stats(file_path: str, block_offsets: List[int]) -> None:
    """
    Print statistics about a BGZF file.

    Args:
        file_path: Path to BGZF file
        block_offsets: List of block offsets from scan_bgzf_blocks()
    """
    file_size = Path(file_path).stat().st_size
    num_blocks = len(block_offsets)
    avg_block_size = file_size / num_blocks if num_blocks > 0 else 0
    estimated_records = estimate_records_from_blocks(num_blocks)

    print(f"BGZF File Statistics:")
    print(f"  File: {file_path}")
    print(f"  Size: {file_size / 1e9:.2f} GB ({file_size:,} bytes)")
    print(f"  Blocks: {num_blocks:,}")
    print(f"  Avg block size: {avg_block_size:.0f} bytes")
    print(f"  Estimated records: ~{estimated_records:,}")
    print(f"  (assuming ~206 records per 64KB block)")


if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:
        print("Usage: python bgzf_scanner.py <bgzf_file> [num_workers]")
        sys.exit(1)

    file_path = sys.argv[1]
    num_workers = int(sys.argv[2]) if len(sys.argv) > 2 else 32

    if not is_bgzf_file(file_path):
        print(f"Error: {file_path} is not a BGZF file")
        sys.exit(1)

    print(f"Scanning BGZF blocks in {file_path}...")
    block_offsets = scan_bgzf_blocks(file_path, progress_interval=1000000)

    print()
    print_bgzf_stats(file_path, block_offsets)

    print()
    print(f"Block distribution for {num_workers} workers:")
    file_size = Path(file_path).stat().st_size
    ranges = get_block_ranges(block_offsets, num_workers, file_size)

    for worker_id, start, end, blocks in ranges[:5]:
        print(f"  Worker {worker_id}: offset {start:,} - {end:,} ({blocks:,} blocks)")
    if len(ranges) > 5:
        print(f"  ... and {len(ranges) - 5} more workers")
