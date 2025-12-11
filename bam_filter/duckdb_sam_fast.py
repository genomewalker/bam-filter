#!/usr/bin/env python3
"""
Ultra-fast SAM.gz to Parquet converter using DuckDB's parallel TSV parser.

Strategy:
1. Use DuckDB's read_csv with null_padding for variable columns
2. Parse optional tags from columns 12+
3. Apply optimized schema (hot/cold tags, packed sequences)
4. Write to Parquet with ZSTD compression

This bypasses HTSlib entirely for maximum speed on SAM.gz files.
"""

import os
import sys
import time
import gzip
from pathlib import Path
from typing import Dict, Any

import duckdb


def count_header_lines(input_file: str) -> int:
    """Count header lines (starting with @) - fast version."""
    count = 0
    open_func = gzip.open if input_file.endswith('.gz') else open

    with open_func(input_file, 'rt') as f:
        for line in f:
            if line.startswith('@'):
                count += 1
            else:
                break

    return count


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


def duckdb_sam_to_parquet_fast(
    input_file: str,
    output_dir: str,
    num_threads: int = 32,
    compression_level: int = 3,
    row_group_size: int = 100000,
) -> Dict[str, Any]:
    """
    Convert SAM.gz to Parquet using DuckDB's parallel TSV parser.

    This is the fastest approach for SAM.gz files because:
    - DuckDB's read_csv is highly optimized and parallel
    - Reads gzip directly without external decompression
    - Writes to Parquet in a single SQL statement

    Args:
        input_file: Input SAM/SAM.gz file
        output_dir: Output directory
        num_threads: Number of DuckDB threads (default: 32)
        compression_level: ZSTD compression level (default: 3)
        row_group_size: Parquet row group size (default: 100000)

    Returns:
        Dictionary with conversion statistics
    """
    print("=" * 80)
    print("DuckDB FAST SAM → Parquet Converter")
    print("=" * 80)
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print(f"Threads: {num_threads}")
    print()

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Count header lines
    print("Counting header lines...")
    header_lines = count_header_lines(input_file)
    print(f"Header lines: {header_lines:,}")

    # Extract header for round-trip capability
    print("Extracting header...")
    extract_header(input_file, output_dir)
    print()

    # Create DuckDB connection
    con = duckdb.connect(":memory:")
    con.execute(f"SET threads TO {num_threads}")
    con.execute("SET preserve_insertion_order = true")

    start_time = time.time()

    # Build the read_csv query
    # Use null_padding to handle variable number of tag columns
    # Use union_by_name to handle columns automatically

    # First, determine column count from a sample
    print("Determining maximum columns...")
    sample_query = f"""
        SELECT max(list_length(string_split(line, '\t'))) as max_cols
        FROM (
            SELECT unnest as line
            FROM (
                SELECT string_split(content, chr(10)) as lines
                FROM read_text('{input_file}', compression='gzip')
            ), unnest(lines)
            WHERE line NOT LIKE '@%' AND line != ''
            LIMIT 10000
        )
    """

    try:
        max_cols = con.execute(sample_query).fetchone()[0]
        print(f"Max columns in sample: {max_cols}")
    except Exception as e:
        print(f"Warning: Could not determine max columns: {e}")
        max_cols = 30  # Reasonable default

    # For SAM format, we need at least 11 columns
    # Columns 12+ are optional tags
    num_cols = max(max_cols, 11)

    # Build column specification - all as VARCHAR for robustness
    columns = {f'column{i}': 'VARCHAR' for i in range(num_cols)}
    columns_str = ', '.join([f"'column{i}': 'VARCHAR'" for i in range(num_cols)])

    print(f"Using {num_cols} columns")
    print()

    # Main conversion query
    # Core SAM fields with proper types + concatenate remaining tags
    output_file = output_path / "alignments.parquet"
    print(f"Converting to {output_file}...")

    # Build tag concatenation expression for columns 11+
    tag_cols = [f"COALESCE(column{i}, '')" for i in range(11, num_cols)]
    if tag_cols:
        tags_expr = "array_to_string(list_filter([" + ", ".join(tag_cols) + "], x -> x != ''), '\t')"
    else:
        tags_expr = "''"

    query = f"""
        COPY (
            SELECT
                row_number() OVER () as read_id,
                column0 AS read_name,
                CAST(column1 AS USMALLINT) AS flag,
                column2 AS ref_name,
                CAST(column3 AS INTEGER) - 1 AS position,  -- Convert to 0-based
                CAST(column4 AS UTINYINT) AS mapq,
                column5 AS cigar,
                column6 AS mate_ref_name,
                CAST(column7 AS INTEGER) - 1 AS mate_position,  -- Convert to 0-based
                CAST(column8 AS INTEGER) AS template_length,
                column9 AS sequence,
                column10 AS quality,
                -- Extract hot tags if present (AS, NM, XS, MD)
                CAST(regexp_extract(tags_raw, 'AS:i:(-?\\d+)', 1) AS INTEGER) AS alignment_score,
                CAST(regexp_extract(tags_raw, 'NM:i:(\\d+)', 1) AS USMALLINT) AS edit_distance,
                CAST(regexp_extract(tags_raw, 'XS:i:(-?\\d+)', 1) AS INTEGER) AS xs_score,
                regexp_extract(tags_raw, 'MD:Z:([^\t]+)', 1) AS md_string,
                regexp_extract(tags_raw, 'RG:Z:([^\t]+)', 1) AS read_group,
                -- Calculate ANI if NM is present
                CASE
                    WHEN regexp_extract(tags_raw, 'NM:i:(\\d+)', 1) IS NOT NULL
                         AND length(column9) > 0
                    THEN (1.0 - CAST(regexp_extract(tags_raw, 'NM:i:(\\d+)', 1) AS FLOAT) / length(column9)) * 100.0
                    ELSE NULL
                END AS ani,
                -- Keep raw tags for lossless storage
                tags_raw
            FROM (
                SELECT
                    *,
                    {tags_expr} AS tags_raw
                FROM read_csv(
                    '{input_file}',
                    delim='\t',
                    header=false,
                    skip={header_lines},
                    compression='gzip',
                    all_varchar=true,
                    ignore_errors=true,
                    null_padding=true,
                    max_line_size=10000000,
                    columns={{{columns_str}}}
                )
            )
            WHERE column0 IS NOT NULL AND column0 != ''
        )
        TO '{output_file}'
        (
            FORMAT PARQUET,
            COMPRESSION 'ZSTD',
            COMPRESSION_LEVEL {compression_level},
            ROW_GROUP_SIZE {row_group_size}
        )
    """

    print("Running conversion...")
    convert_start = time.time()
    con.execute(query)
    convert_time = time.time() - convert_start

    total_time = time.time() - start_time

    # Get output stats
    output_size = output_file.stat().st_size

    # Count records
    print("Counting output records...")
    count_result = con.execute(f"SELECT COUNT(*) FROM '{output_file}'").fetchone()[0]

    print()
    print("=" * 80)
    print("DUCKDB FAST CONVERSION COMPLETE")
    print("=" * 80)
    print(f"Total records: {count_result:,}")
    print(f"Conversion time: {convert_time:.1f}s")
    print(f"Total time: {total_time:.1f}s ({total_time/60:.2f} minutes)")
    print(f"Throughput: {count_result/total_time:,.0f} records/sec")
    print(f"Throughput: {count_result/total_time/1e6:.3f} M records/sec")
    print(f"Output size: {output_size/1e9:.2f} GB ({output_size/count_result:.1f} bytes/record)")
    print()

    return {
        'total_records': count_result,
        'convert_time_seconds': convert_time,
        'total_time_seconds': total_time,
        'throughput_records_per_sec': count_result / total_time if total_time > 0 else 0,
        'output_size_bytes': output_size,
    }


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="DuckDB Fast SAM to Parquet Converter")
    parser.add_argument('input_file', help="Input SAM/SAM.gz file")
    parser.add_argument('output_dir', help="Output directory")
    parser.add_argument('-t', '--threads', type=int, default=32, help="Number of threads")
    parser.add_argument('--compression-level', type=int, default=3, help="ZSTD compression level")
    parser.add_argument('--row-group-size', type=int, default=100000, help="Parquet row group size")

    args = parser.parse_args()

    stats = duckdb_sam_to_parquet_fast(
        args.input_file,
        args.output_dir,
        num_threads=args.threads,
        compression_level=args.compression_level,
        row_group_size=args.row_group_size,
    )
