#!/usr/bin/env python3
"""
Prototype: High-performance SAM.gz reader using DuckDB's parallel CSV reader.

SAM format: TAB-separated with 11 mandatory columns + optional tags
@HD, @SQ, @RG lines need to be skipped (header)

Columns:
1. QNAME - Query template NAME
2. FLAG  - bitwise FLAG
3. RNAME - Reference sequence NAME
4. POS   - 1-based leftmost mapping POSition
5. MAPQ  - MAPping Quality
6. CIGAR - CIGAR string
7. RNEXT - Ref. name of the mate/next read
8. PNEXT - Position of the mate/next read
9. TLEN  - observed Template LENgth
10. SEQ   - segment SEQuence
11. QUAL  - ASCII of Phred-scaled base QUALity+33
12+ TAGS - Optional fields (AS:i:, NM:i:, etc.)

DuckDB can parallelize CSV reading, but sam.gz limits parallelism due to
gzip's sequential nature. BGZF helps but DuckDB doesn't natively seek BGZF.

Approach: Decompress BGZF blocks in parallel, then use DuckDB on decompressed text.
"""

import sys
import time
import subprocess
import tempfile
from pathlib import Path

# Check if DuckDB available
try:
    import duckdb
    HAS_DUCKDB = True
except ImportError:
    HAS_DUCKDB = False


def count_header_lines(sam_gz_path: str) -> int:
    """Count lines starting with @ in a sam.gz file."""
    result = subprocess.run(
        f"zcat {sam_gz_path} | head -10000 | grep -c '^@' || true",
        shell=True, capture_output=True, text=True
    )
    return int(result.stdout.strip() or 0)


def test_duckdb_sam_reader(sam_gz_path: str, max_records: int = 1000000):
    """
    Test DuckDB's CSV reader on SAM.gz file.

    Note: DuckDB reads gzipped files sequentially (single-threaded decompression).
    For true parallelism, we'd need to decompress BGZF blocks first.
    """
    if not HAS_DUCKDB:
        print("DuckDB not installed. Install with: pip install duckdb")
        return

    print(f"Testing DuckDB SAM reader on: {sam_gz_path}")

    # Count header lines to skip
    print("Counting header lines...")
    header_lines = count_header_lines(sam_gz_path)
    print(f"  Header lines: {header_lines}")

    # Define SAM schema (11 mandatory columns + tags as text)
    sam_columns = {
        'qname': 'VARCHAR',
        'flag': 'USMALLINT',  # uint16
        'rname': 'VARCHAR',
        'pos': 'INTEGER',
        'mapq': 'UTINYINT',  # uint8
        'cigar': 'VARCHAR',
        'rnext': 'VARCHAR',
        'pnext': 'INTEGER',
        'tlen': 'INTEGER',
        'seq': 'VARCHAR',
        'qual': 'VARCHAR',
    }

    col_names = list(sam_columns.keys())
    col_types = list(sam_columns.values())

    # Connect to DuckDB (in-memory)
    con = duckdb.connect(':memory:')

    print(f"\nReading with DuckDB (skip_rows={header_lines})...")
    start_time = time.time()

    # DuckDB read_csv with SAM settings
    # Note: DuckDB will be slower on gzip because it can't parallelize decompression
    query = f"""
    SELECT COUNT(*) as cnt
    FROM read_csv(
        '{sam_gz_path}',
        delim='\t',
        header=false,
        skip={header_lines},
        columns={{
            'qname': 'VARCHAR',
            'flag': 'USMALLINT',
            'rname': 'VARCHAR',
            'pos': 'INTEGER',
            'mapq': 'UTINYINT',
            'cigar': 'VARCHAR',
            'rnext': 'VARCHAR',
            'pnext': 'INTEGER',
            'tlen': 'INTEGER',
            'seq': 'VARCHAR',
            'qual': 'VARCHAR'
        }},
        ignore_errors=true
    )
    WHERE pos > 0
    """

    try:
        result = con.execute(query).fetchone()
        count = result[0] if result else 0
        elapsed = time.time() - start_time

        print(f"\nResults:")
        print(f"  Records: {count:,}")
        print(f"  Time: {elapsed:.1f}s")
        print(f"  Throughput: {count/elapsed:,.0f} rec/s")

        # Sample some data
        sample_query = f"""
        SELECT qname, flag, rname, pos, mapq, cigar, tlen
        FROM read_csv(
            '{sam_gz_path}',
            delim='\t',
            header=false,
            skip={header_lines},
            columns={{
                'qname': 'VARCHAR',
                'flag': 'USMALLINT',
                'rname': 'VARCHAR',
                'pos': 'INTEGER',
                'mapq': 'UTINYINT',
                'cigar': 'VARCHAR',
                'rnext': 'VARCHAR',
                'pnext': 'INTEGER',
                'tlen': 'INTEGER',
                'seq': 'VARCHAR',
                'qual': 'VARCHAR'
            }},
            ignore_errors=true
        )
        WHERE pos > 0
        LIMIT 5
        """

        print("\nSample records:")
        sample = con.execute(sample_query).fetchdf()
        print(sample.to_string())

    except Exception as e:
        print(f"Error: {e}")

    finally:
        con.close()


def test_duckdb_decompressed(sam_gz_path: str):
    """
    Test: Decompress to temporary file, then use DuckDB's parallel reader.

    This should be faster because DuckDB can parallelize reading uncompressed data.
    """
    if not HAS_DUCKDB:
        print("DuckDB not installed")
        return

    print(f"\nTesting DuckDB with decompression to temp file...")
    print("(This allows DuckDB to parallelize reading)")

    # Count header lines
    header_lines = count_header_lines(sam_gz_path)

    # Create temp file
    with tempfile.NamedTemporaryFile(suffix='.sam', delete=False) as tmp:
        tmp_path = tmp.name

    try:
        # Decompress with pigz for parallel decompression
        print(f"Decompressing to {tmp_path}...")
        start_decomp = time.time()
        subprocess.run(
            f"zcat {sam_gz_path} > {tmp_path}",
            shell=True, check=True
        )
        decomp_time = time.time() - start_decomp
        print(f"  Decompression time: {decomp_time:.1f}s")

        # Get file size
        size_mb = Path(tmp_path).stat().st_size / 1e6
        print(f"  Uncompressed size: {size_mb:.1f} MB")

        # Now read with DuckDB (parallel)
        con = duckdb.connect(':memory:')
        con.execute("SET threads=8")

        print("\nReading with DuckDB (parallel)...")
        start_time = time.time()

        query = f"""
        SELECT COUNT(*) as cnt
        FROM read_csv(
            '{tmp_path}',
            delim='\t',
            header=false,
            skip={header_lines},
            columns={{
                'qname': 'VARCHAR',
                'flag': 'USMALLINT',
                'rname': 'VARCHAR',
                'pos': 'INTEGER',
                'mapq': 'UTINYINT',
                'cigar': 'VARCHAR',
                'rnext': 'VARCHAR',
                'pnext': 'INTEGER',
                'tlen': 'INTEGER',
                'seq': 'VARCHAR',
                'qual': 'VARCHAR'
            }},
            parallel=true,
            ignore_errors=true
        )
        WHERE pos > 0
        """

        result = con.execute(query).fetchone()
        count = result[0] if result else 0
        read_time = time.time() - start_time
        total_time = decomp_time + read_time

        print(f"\nResults:")
        print(f"  Records: {count:,}")
        print(f"  Read time: {read_time:.1f}s")
        print(f"  Total time: {total_time:.1f}s")
        print(f"  Read throughput: {count/read_time:,.0f} rec/s")
        print(f"  Total throughput: {count/total_time:,.0f} rec/s")

        con.close()

    finally:
        Path(tmp_path).unlink(missing_ok=True)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python test_duckdb_sam_reader.py <input.sam.gz>")
        sys.exit(1)

    sam_gz = sys.argv[1]

    # Test 1: DuckDB on gzipped file (sequential decompression)
    test_duckdb_sam_reader(sam_gz)

    # Test 2: Decompress first, then parallel DuckDB
    print("\n" + "="*60)
    test_duckdb_decompressed(sam_gz)
