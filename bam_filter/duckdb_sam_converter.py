#!/usr/bin/env python3
"""
DuckDB-based SAM/BAM to Parquet Converter.

Strategy: Use DuckDB's highly optimized parallel TSV parser to read SAM files
directly, then write to Parquet. This avoids the overhead of HTSlib's
sequential SAM parsing.

SAM format:
- Lines starting with @ are header lines
- Data lines have 11+ tab-separated columns:
  1. QNAME  - query name
  2. FLAG   - bitwise flag
  3. RNAME  - reference name
  4. POS    - 1-based position
  5. MAPQ   - mapping quality
  6. CIGAR  - CIGAR string
  7. RNEXT  - mate reference name
  8. PNEXT  - mate position
  9. TLEN   - template length
  10. SEQ   - sequence
  11. QUAL  - quality string
  12+ - Optional tags in TAG:TYPE:VALUE format
"""

import os
import sys
import time
import gzip
from pathlib import Path
from typing import Dict, Any, Optional, List

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq


# SAM mandatory columns
SAM_COLUMNS = [
    ('QNAME', 'VARCHAR'),
    ('FLAG', 'INTEGER'),
    ('RNAME', 'VARCHAR'),
    ('POS', 'INTEGER'),
    ('MAPQ', 'INTEGER'),
    ('CIGAR', 'VARCHAR'),
    ('RNEXT', 'VARCHAR'),
    ('PNEXT', 'INTEGER'),
    ('TLEN', 'INTEGER'),
    ('SEQ', 'VARCHAR'),
    ('QUAL', 'VARCHAR'),
]

# Common optional tags we want to extract
OPTIONAL_TAGS = ['NM', 'AS', 'XS', 'MD', 'NH', 'RG', 'BC', 'CB', 'UB']


def extract_header(input_file: str, output_dir: str) -> str:
    """Extract SAM header to a separate file."""
    header_file = Path(output_dir) / "header.sam"

    open_func = gzip.open if input_file.endswith('.gz') else open

    with open_func(input_file, 'rt') as f_in:
        with open(header_file, 'w') as f_out:
            for line in f_in:
                if line.startswith('@'):
                    f_out.write(line)
                else:
                    break

    return str(header_file)


def count_header_lines(input_file: str) -> int:
    """Count header lines (starting with @)."""
    count = 0
    open_func = gzip.open if input_file.endswith('.gz') else open

    with open_func(input_file, 'rt') as f:
        for line in f:
            if line.startswith('@'):
                count += 1
            else:
                break

    return count


def create_sam_schema() -> pa.Schema:
    """Create PyArrow schema for SAM alignment records."""
    return pa.schema([
        pa.field('qname', pa.string()),
        pa.field('flag', pa.uint16()),
        pa.field('reference', pa.string()),
        pa.field('pos', pa.int32()),
        pa.field('mapq', pa.uint8()),
        pa.field('cigar', pa.string()),
        pa.field('rnext', pa.string()),
        pa.field('pnext', pa.int32()),
        pa.field('tlen', pa.int32()),
        pa.field('seq', pa.string()),
        pa.field('qual', pa.string()),
        # Common optional tags
        pa.field('NM', pa.int32()),  # edit distance
        pa.field('AS', pa.int32()),  # alignment score
        pa.field('XS', pa.int32()),  # suboptimal score
        pa.field('MD', pa.string()), # mismatching positions
        pa.field('RG', pa.string()), # read group
    ])


def duckdb_convert_sam_to_parquet(
    input_file: str,
    output_dir: str,
    batch_size: int = 1000000,
    num_threads: int = 32,
    compression_level: int = 3,
) -> Dict[str, Any]:
    """
    Convert SAM to Parquet using DuckDB's parallel TSV parser.

    Args:
        input_file: Input SAM/SAM.gz file
        output_dir: Output directory
        batch_size: Records per batch
        num_threads: Number of DuckDB threads
        compression_level: ZSTD compression level

    Returns:
        Dictionary with conversion statistics
    """
    print("=" * 80)
    print("DuckDB-based SAM → Parquet Converter")
    print("=" * 80)
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print(f"Threads: {num_threads}")
    print()

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Count header lines to skip
    print("Counting header lines...")
    header_lines = count_header_lines(input_file)
    print(f"Header lines: {header_lines:,}")

    # Extract header
    print("Extracting header...")
    extract_header(input_file, output_dir)
    print()

    # Create DuckDB connection
    con = duckdb.connect(":memory:")
    con.execute(f"SET threads TO {num_threads}")

    # DuckDB can read gzip directly and skip header lines
    # We read all columns as VARCHAR first, then cast
    print("Reading SAM file with DuckDB...")
    start_time = time.time()

    # Create a view that reads the SAM file
    # Note: SAM files have variable columns due to optional tags
    # We read first 11 columns (mandatory) plus concatenate remaining as tags
    query = f"""
        SELECT
            column0 AS qname,
            CAST(column1 AS INTEGER) AS flag,
            column2 AS reference,
            CAST(column3 AS INTEGER) AS pos,
            CAST(column4 AS INTEGER) AS mapq,
            column5 AS cigar,
            column6 AS rnext,
            CAST(column7 AS INTEGER) AS pnext,
            CAST(column8 AS INTEGER) AS tlen,
            column9 AS seq,
            column10 AS qual
        FROM read_csv(
            '{input_file}',
            delim='\t',
            header=false,
            skip={header_lines},
            compression='gzip',
            all_varchar=true,
            ignore_errors=true,
            max_line_size=10000000,
            columns={{
                'column0': 'VARCHAR',
                'column1': 'VARCHAR',
                'column2': 'VARCHAR',
                'column3': 'VARCHAR',
                'column4': 'VARCHAR',
                'column5': 'VARCHAR',
                'column6': 'VARCHAR',
                'column7': 'VARCHAR',
                'column8': 'VARCHAR',
                'column9': 'VARCHAR',
                'column10': 'VARCHAR'
            }}
        )
    """

    # First, count total records
    print("Counting total records...")
    count_query = f"""
        SELECT COUNT(*) FROM read_csv(
            '{input_file}',
            delim='\t',
            header=false,
            skip={header_lines},
            compression='gzip',
            all_varchar=true,
            ignore_errors=true,
            max_line_size=10000000
        )
    """

    count_start = time.time()
    total_records = con.execute(count_query).fetchone()[0]
    count_time = time.time() - count_start
    print(f"Total records: {total_records:,}")
    print(f"Count time: {count_time:.1f}s ({total_records/count_time:,.0f} rec/sec)")
    print()

    # Write directly to parquet using DuckDB
    output_file = output_path / "alignments.parquet"
    print(f"Writing to {output_file}...")

    write_query = f"""
        COPY ({query})
        TO '{output_file}'
        (FORMAT PARQUET, COMPRESSION 'ZSTD', COMPRESSION_LEVEL {compression_level})
    """

    write_start = time.time()
    con.execute(write_query)
    write_time = time.time() - write_start

    total_time = time.time() - start_time

    # Get output file size
    output_size = output_file.stat().st_size

    print()
    print("=" * 80)
    print("DUCKDB CONVERSION COMPLETE")
    print("=" * 80)
    print(f"Total records: {total_records:,}")
    print(f"Count time: {count_time:.1f}s")
    print(f"Write time: {write_time:.1f}s")
    print(f"Total time: {total_time:.1f}s ({total_time/60:.2f} minutes)")
    print(f"Throughput: {total_records/total_time:,.0f} records/sec")
    print(f"Throughput: {total_records/total_time/1e6:.3f} M records/sec")
    print(f"Output size: {output_size/1e9:.2f} GB")
    print()

    return {
        'total_records': total_records,
        'count_time_seconds': count_time,
        'write_time_seconds': write_time,
        'total_time_seconds': total_time,
        'throughput_records_per_sec': total_records / total_time if total_time > 0 else 0,
        'output_size_bytes': output_size,
    }


def duckdb_convert_sam_to_parquet_batched(
    input_file: str,
    output_dir: str,
    batch_size: int = 10000000,
    num_threads: int = 32,
    compression_level: int = 3,
) -> Dict[str, Any]:
    """
    Convert SAM to Parquet using DuckDB with batched processing.

    This version processes the file in batches, writing multiple smaller
    parquet files that can be read together.

    Args:
        input_file: Input SAM/SAM.gz file
        output_dir: Output directory
        batch_size: Records per parquet file
        num_threads: Number of DuckDB threads
        compression_level: ZSTD compression level

    Returns:
        Dictionary with conversion statistics
    """
    print("=" * 80)
    print("DuckDB-based SAM → Parquet Converter (Batched)")
    print("=" * 80)
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print(f"Threads: {num_threads}")
    print(f"Batch size: {batch_size:,}")
    print()

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Count header lines to skip
    print("Counting header lines...")
    header_lines = count_header_lines(input_file)
    print(f"Header lines: {header_lines:,}")

    # Extract header
    print("Extracting header...")
    extract_header(input_file, output_dir)
    print()

    # Create DuckDB connection
    con = duckdb.connect(":memory:")
    con.execute(f"SET threads TO {num_threads}")
    con.execute("SET preserve_insertion_order = true")

    start_time = time.time()

    # Use DuckDB's EXPORT feature with row_group_size for batching
    # This writes a single parquet file with multiple row groups
    query = f"""
        SELECT
            column0 AS qname,
            CAST(column1 AS USMALLINT) AS flag,
            column2 AS reference,
            CAST(column3 AS INTEGER) AS pos,
            CAST(column4 AS UTINYINT) AS mapq,
            column5 AS cigar,
            column6 AS rnext,
            CAST(column7 AS INTEGER) AS pnext,
            CAST(column8 AS INTEGER) AS tlen,
            column9 AS seq,
            column10 AS qual
        FROM read_csv(
            '{input_file}',
            delim='\t',
            header=false,
            skip={header_lines},
            compression='gzip',
            all_varchar=true,
            ignore_errors=true,
            max_line_size=10000000,
            columns={{
                'column0': 'VARCHAR',
                'column1': 'VARCHAR',
                'column2': 'VARCHAR',
                'column3': 'VARCHAR',
                'column4': 'VARCHAR',
                'column5': 'VARCHAR',
                'column6': 'VARCHAR',
                'column7': 'VARCHAR',
                'column8': 'VARCHAR',
                'column9': 'VARCHAR',
                'column10': 'VARCHAR'
            }}
        )
    """

    # Write with partitioning by reference for efficient queries
    # Or use row_group_size for batching within single file
    output_file = output_path / "alignments.parquet"
    print(f"Converting to {output_file}...")

    write_query = f"""
        COPY ({query})
        TO '{output_file}'
        (
            FORMAT PARQUET,
            COMPRESSION 'ZSTD',
            COMPRESSION_LEVEL {compression_level},
            ROW_GROUP_SIZE {batch_size}
        )
    """

    con.execute(write_query)

    total_time = time.time() - start_time

    # Get stats
    output_size = output_file.stat().st_size

    # Count records in output
    count_result = con.execute(f"SELECT COUNT(*) FROM '{output_file}'").fetchone()[0]

    print()
    print("=" * 80)
    print("DUCKDB BATCHED CONVERSION COMPLETE")
    print("=" * 80)
    print(f"Total records: {count_result:,}")
    print(f"Total time: {total_time:.1f}s ({total_time/60:.2f} minutes)")
    print(f"Throughput: {count_result/total_time:,.0f} records/sec")
    print(f"Throughput: {count_result/total_time/1e6:.3f} M records/sec")
    print(f"Output size: {output_size/1e9:.2f} GB")
    print()

    return {
        'total_records': count_result,
        'total_time_seconds': total_time,
        'throughput_records_per_sec': count_result / total_time if total_time > 0 else 0,
        'output_size_bytes': output_size,
    }


def duckdb_convert_sam_partitioned(
    input_file: str,
    output_dir: str,
    num_threads: int = 32,
    compression_level: int = 3,
) -> Dict[str, Any]:
    """
    Convert SAM to Parquet with partitioning by reference (chromosome).

    This creates a directory structure:
    output_dir/
        reference=chr1/part_0.parquet
        reference=chr2/part_0.parquet
        ...

    This is optimal for queries filtering by reference.
    """
    print("=" * 80)
    print("DuckDB-based SAM → Parquet Converter (Partitioned by Reference)")
    print("=" * 80)
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print(f"Threads: {num_threads}")
    print()

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    parquet_dir = output_path / "by_reference"
    parquet_dir.mkdir(exist_ok=True)

    # Count header lines
    print("Counting header lines...")
    header_lines = count_header_lines(input_file)
    print(f"Header lines: {header_lines:,}")

    # Extract header
    print("Extracting header...")
    extract_header(input_file, output_dir)
    print()

    # Create DuckDB connection
    con = duckdb.connect(":memory:")
    con.execute(f"SET threads TO {num_threads}")

    start_time = time.time()

    query = f"""
        SELECT
            column0 AS qname,
            CAST(column1 AS USMALLINT) AS flag,
            column2 AS reference,
            CAST(column3 AS INTEGER) AS pos,
            CAST(column4 AS UTINYINT) AS mapq,
            column5 AS cigar,
            column6 AS rnext,
            CAST(column7 AS INTEGER) AS pnext,
            CAST(column8 AS INTEGER) AS tlen,
            column9 AS seq,
            column10 AS qual
        FROM read_csv(
            '{input_file}',
            delim='\t',
            header=false,
            skip={header_lines},
            compression='gzip',
            all_varchar=true,
            ignore_errors=true,
            max_line_size=10000000,
            columns={{
                'column0': 'VARCHAR',
                'column1': 'VARCHAR',
                'column2': 'VARCHAR',
                'column3': 'VARCHAR',
                'column4': 'VARCHAR',
                'column5': 'VARCHAR',
                'column6': 'VARCHAR',
                'column7': 'VARCHAR',
                'column8': 'VARCHAR',
                'column9': 'VARCHAR',
                'column10': 'VARCHAR'
            }}
        )
    """

    print("Converting with Hive partitioning by reference...")

    write_query = f"""
        COPY ({query})
        TO '{parquet_dir}'
        (
            FORMAT PARQUET,
            COMPRESSION 'ZSTD',
            COMPRESSION_LEVEL {compression_level},
            PARTITION_BY (reference),
            OVERWRITE_OR_IGNORE true
        )
    """

    con.execute(write_query)

    total_time = time.time() - start_time

    # Count output files and total size
    parquet_files = list(parquet_dir.rglob("*.parquet"))
    total_size = sum(f.stat().st_size for f in parquet_files)

    # Count records
    count_result = con.execute(f"SELECT COUNT(*) FROM '{parquet_dir}/**/*.parquet'").fetchone()[0]

    print()
    print("=" * 80)
    print("DUCKDB PARTITIONED CONVERSION COMPLETE")
    print("=" * 80)
    print(f"Total records: {count_result:,}")
    print(f"Partitions: {len(parquet_files)}")
    print(f"Total time: {total_time:.1f}s ({total_time/60:.2f} minutes)")
    print(f"Throughput: {count_result/total_time:,.0f} records/sec")
    print(f"Throughput: {count_result/total_time/1e6:.3f} M records/sec")
    print(f"Output size: {total_size/1e9:.2f} GB")
    print()

    return {
        'total_records': count_result,
        'num_partitions': len(parquet_files),
        'total_time_seconds': total_time,
        'throughput_records_per_sec': count_result / total_time if total_time > 0 else 0,
        'output_size_bytes': total_size,
    }


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="DuckDB-based SAM to Parquet Converter")
    parser.add_argument('input_file', help="Input SAM/SAM.gz file")
    parser.add_argument('output_dir', help="Output directory")
    parser.add_argument('-t', '--threads', type=int, default=32, help="Number of threads")
    parser.add_argument('--mode', choices=['simple', 'batched', 'partitioned'],
                        default='batched', help="Conversion mode")
    parser.add_argument('--batch-size', type=int, default=10000000,
                        help="Batch size for batched mode")
    parser.add_argument('--compression-level', type=int, default=3,
                        help="ZSTD compression level")

    args = parser.parse_args()

    if args.mode == 'simple':
        stats = duckdb_convert_sam_to_parquet(
            args.input_file,
            args.output_dir,
            num_threads=args.threads,
            compression_level=args.compression_level,
        )
    elif args.mode == 'batched':
        stats = duckdb_convert_sam_to_parquet_batched(
            args.input_file,
            args.output_dir,
            batch_size=args.batch_size,
            num_threads=args.threads,
            compression_level=args.compression_level,
        )
    else:  # partitioned
        stats = duckdb_convert_sam_partitioned(
            args.input_file,
            args.output_dir,
            num_threads=args.threads,
            compression_level=args.compression_level,
        )
