import duckdb
import logging
import os
import psutil
from pathlib import Path
import json
import numpy as np
from numba import njit
import shutil  # Added for shutil.which and disk_usage
import subprocess  # Added for subprocess.Popen
import tempfile  # Added for temporary named pipe creation
import bz2
import sys  # Import sys
from bam_filter import __version__  # Import version
from .parquet_utils import calculate_optimal_row_group_size, get_column_info, ColumnInfo


log = logging.getLogger("my_logger")

# --- DNA Reverse Complement Logic ---
complement_map = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def dna_revcomp(seq):
    """Reverse complement a DNA sequence."""
    if seq is None:
        return None
    return seq.translate(complement_map)[::-1]


def str_reverse(s):
    """Reverse a string."""
    if s is None:
        return None
    return s[::-1]


# --- End DNA Reverse Complement Logic ---


# --- Moved get_file_extension function here ---
def get_file_extension(file_path):
    """
    Determines the base extension and compression type of a file based on its suffix(es).
    Handles extensions like .bam, .sam, .sam.gz, .sam.bz2, .parquet.

    Args:
        file_path (str): The path to the file.

    Returns:
        tuple: (base_extension, compression_type)
               e.g., ('.bam', None), ('.sam', 'gz'), ('.sam', 'bz2'), ('.sam', None), ('.parquet', None)
               Returns (last_suffix, None) or (None, None) if the extension pattern is unrecognized.
    """
    path = Path(file_path)
    name = path.name.lower()  # Use lower case for consistent matching

    # Check common patterns first
    if name.endswith(".bam"):
        return ".bam", None
    elif name.endswith(".sam.gz"):
        return ".sam", "gz"
    elif name.endswith(".sam.bz2"):
        return ".sam", "bz2"
    elif name.endswith(".sam"):
        return ".sam", None
    elif name.endswith(".parquet"):
        # This identifies a file ending in .parquet.
        # Directory detection is handled in convert.py
        return ".parquet", None
    else:
        # Fallback for less common or unknown extensions
        suffixes = path.suffixes
        if len(suffixes) >= 2:
            base_ext = suffixes[-2].lower()
            compression_ext = suffixes[-1].lower()
            if compression_ext == ".gz":
                return base_ext, "gz"
            elif compression_ext == ".bz2":
                return base_ext, "bz2"
            else:
                # If the last suffix isn't compression, return the full suffix as base
                return "".join(suffixes).lower(), None
        elif len(suffixes) == 1:
            # Single suffix, assume no compression
            return suffixes[0].lower(), None
        else:
            # No suffix found
            return None, None


# --- End moved function ---


# Define ONLY the tags essential for score/ANI calculation to extract into columns
COMMON_TAGS = {
    "ZA": "FLOAT",  # ANI value (Average Nucleotide Identity) - Custom tag
    "ZS": "FLOAT",  # Shifted score - Custom tag
}


def export_optimized_table(con, table_name, output_dir_str, compression_level=11):
    """
    Export a table to Parquet format with optimized settings.
    For the 'alignments' table, Hive partitioning by 'partition_id' is applied.
    Sorting is skipped to improve export performance.

    Args:
        con: DuckDB connection object.
        table_name (str): Name of the table to export.
        output_dir_str (str): Path to the output base directory.
        compression_level (int): ZSTD compression level (1-22).

    Returns:
        bool: True if export was successful, False otherwise.
    """
    output_dir = Path(output_dir_str)
    partition_options = ""

    # Determine the export path and partitioning options
    if table_name == "alignments":
        # For alignments table, output to a subdirectory named 'alignments'
        # This path will become the base for Hive partitioning.
        export_path = output_dir / table_name
        # Check if partition_id column exists before attempting to partition
        try:
            col_info = con.execute(f"DESCRIBE {table_name}").fetchall()
            if any(col[0] == "partition_id" for col in col_info):
                partition_options = ", PARTITION_BY (partition_id)"
                log.info(
                    f"Hive partitioning by 'partition_id' will be applied to '{table_name}'. Output will be a directory: {export_path}"
                )
            else:
                log.warning(
                    f"'partition_id' column not found in '{table_name}'. Exporting without Hive partitioning to a single file."
                )
                export_path = output_dir / f"{table_name}.parquet"
        except Exception as e:
            log.warning(
                f"Could not describe table {table_name} to check for partition_id: {e}. Exporting without Hive partitioning to a single file."
            )
            export_path = output_dir / f"{table_name}.parquet"
    else:
        export_path = output_dir / f"{table_name}.parquet"

    try:
        # Ensure the output directory exists
        if table_name != "alignments" or not partition_options:
            export_path.parent.mkdir(parents=True, exist_ok=True)
        else:
            output_dir.mkdir(parents=True, exist_ok=True)

        # OPTIMIZATION: Calculate appropriate row group size
        try:
            table_schema = con.execute(f"DESCRIBE {table_name}").fetchall()
            columns_info_dict = {}
            for col in table_schema:
                col_name = col[0]
                col_type = col[1]
                columns_info_dict[col_name] = get_column_info(col_type)

            row_count = con.execute(f'SELECT COUNT(*) FROM "{table_name}"').fetchone()[
                0
            ]
            
            # Calculate optimal row group size based on data characteristics
            row_group_size = calculate_optimal_row_group_size(
                con, table_name, columns_info_dict
            )
            
            # Fix: Cap the row group size to ensure we have multiple row groups for better query performance
            # For small tables, use the row count itself, for larger tables aim for 4-8 row groups
            if row_count > 0:
                if row_count <= 100000:  # Very small table
                    # Just use the row count as is for small tables
                    adjusted_row_group_size = row_count
                else:
                    # Calculate target row groups - bigger tables need more row groups
                    target_row_groups = max(4, min(8, row_count // 2000000 + 2))
                    adjusted_row_group_size = min(
                        row_group_size,  # The calculated optimal size
                        row_count,       # Never exceed total rows
                        (row_count + target_row_groups - 1) // target_row_groups  # Aim for target number of groups
                    )
                
                # Provide detailed debugging information
                log.debug(f"""
                Row group size calculation for table '{table_name}':
                - Row size: {columns_info_dict['__avg_row_size__'] if '__avg_row_size__' in columns_info_dict else 'Unknown'} bytes
                - Total rows: {row_count:,}
                - Target group size: {row_group_size:,} (uncapped)
                - Adjusted rows per group: {adjusted_row_group_size:,}
                - Target row groups: {(row_count + adjusted_row_group_size - 1) // adjusted_row_group_size if adjusted_row_group_size > 0 else 0}
                """)
                
                row_group_size = adjusted_row_group_size
            else:
                # For empty tables, use a default
                row_group_size = 100000
                log.warning(f"Table '{table_name}' appears to be empty. Using default row group size.")
                
            row_group_size_option = f", ROW_GROUP_SIZE {row_group_size}"
            log.info(
                f"Using optimal calculated row group size of {row_group_size:,} rows for table '{table_name}' "
                f"({(row_count + row_group_size - 1) // row_group_size if row_group_size > 0 else 0} row groups)"
            )
        except Exception as e:
            log.warning(
                f"Could not determine optimal row group size: {e}. Using default."
            )
            row_group_size_option = ""

        # Construct the SELECT statement without ORDER BY for better performance
        # Special case for tags table - exclude tag_hash and tag_bitmap columns
        if table_name == "tags":
            select_statement = 'SELECT tag_id, raw_tags FROM "tags"'
            log.info(f"Exporting tags table without tag_hash and tag_bitmap columns")
        else:
            select_statement = f'SELECT * FROM "{table_name}"'

        # Skip sorting for faster export performance
        if table_name == "alignments":
            log.info(f"Skipping sorting for '{table_name}' table to improve export speed")
        elif table_name in {"refs", "reads", "cigars", "tags"}:
            log.info(f"Skipping sorting for '{table_name}' table to improve export speed")

        # Create export query with ZSTD compression
        copy_query = f"""
        COPY ({select_statement}) 
        TO '{str(export_path)}' 
        (
            FORMAT PARQUET, 
            PARQUET_VERSION v2,
            COMPRESSION_LEVEL {compression_level},
            COMPRESSION ZSTD{row_group_size_option}{partition_options},
            OVERWRITE_OR_IGNORE TRUE
        );
        """

        log.info(
            f"Exporting table '{table_name}' to '{export_path}' with ZSTD compression..."
        )
        con.execute(copy_query)
        log.info(f"Successfully exported table '{table_name}' to '{export_path}'.")
        return True
    except Exception as e:
        log.error(f"Failed to export table '{table_name}' to '{export_path}': {e}")
        log.debug(f"Query used: {copy_query if 'copy_query' in locals() else 'N/A'}")
        return False


def create_tag_extract_sql(source_column_name="_full_raw_tags"):
    """Create SQL expressions to extract tag values from a source column."""
    tag_extracts = []
    for tag, tag_type in COMMON_TAGS.items():
        col_name = f"tag_{tag}"
        # Create case statement to handle missing tags
        tag_extracts.append(
            f"""
            CASE 
                WHEN regexp_matches({source_column_name}, '{tag}:[AifZHB]:[^\\t]+') THEN
                    REGEXP_EXTRACT({source_column_name}, '{tag}:[AifZHB]:([^\\t]+)', 1)
                ELSE NULL
            END as {col_name}
        """
        )
    return ",\n".join(tag_extracts)


def extract_sam_headers(file_path, threads=1):
    """
    Extract SAM headers from a file using samtools view -H.
    Works with both BAM and SAM files (compressed or uncompressed).

    Args:
        file_path (str): Path to the SAM/BAM file
        threads (int): Number of threads to use for samtools

    Returns:
        list: List of header lines
    """
    cmd = ["samtools", "view", "-H"]

    # Add threading if more than 1 thread requested
    if threads > 1:
        cmd.extend(["-@", str(max(1, threads - 1))])

    cmd.append(file_path)

    log.info(f"Extracting headers using: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        header_lines = result.stdout.splitlines()
        log.info(f"Extracted {len(header_lines)} header lines from {file_path}")
        return header_lines
    except subprocess.CalledProcessError as e:
        log.error(f"Failed to extract headers with samtools: {e}")
        log.error(f"stderr: {e.stderr}")
        return []
    except Exception as e:
        log.error(f"Unexpected error extracting headers: {e}")
        return []


def create_header_tables_from_raw_headers(con, header_lines):
    """
    Create header table from raw header lines extracted with samtools.
    Uses named pipe in /dev/shm (if available) for optimal performance.

    Args:
        con: DuckDB connection
        header_lines (list): List of header lines

    Returns:
        int: Number of header lines processed
    """
    if not header_lines:
        log.warning("No header lines provided to process")
        return 0

    # Create a named pipe in /dev/shm if available, otherwise in tmp
    import os
    import uuid

    pipe_path = None
    use_pipe = False

    try:
        # Check if /dev/shm is available
        if os.path.exists("/dev/shm") and os.access("/dev/shm", os.W_OK):
            pipe_path = f"/dev/shm/bam_filter_header_{uuid.uuid4()}"
            use_pipe = True
        else:
            # Fall back to regular tmp dir
            pipe_path = os.path.join(
                tempfile.gettempdir(), f"bam_filter_header_{uuid.uuid4()}"
            )
            use_pipe = True

        if use_pipe:
            # Create named pipe
            os.mkfifo(pipe_path)
            log.debug(f"Created named pipe for headers at {pipe_path}")

            # Start a background thread to write header lines to the pipe
            import threading

            def write_headers():
                with open(pipe_path, "w") as f:
                    for line in header_lines:
                        f.write(f"{line}\n")

            thread = threading.Thread(target=write_headers)
            thread.daemon = True
            thread.start()

            # Create the header table directly from the named pipe
            query = f"""
            CREATE TABLE header AS
            SELECT 
                header_line,
                -- More efficient pattern matching for header types
                -- Extract once and compare instead of multiple substrings
                CASE
                    WHEN STARTS_WITH(SUBSTR(header_line, 2, 2), 'HD') THEN 'HD'
                    WHEN STARTS_WITH(SUBSTR(header_line, 2, 2), 'SQ') THEN 'SQ'
                    WHEN STARTS_WITH(SUBSTR(header_line, 2, 2), 'RG') THEN 'RG'
                    WHEN STARTS_WITH(SUBSTR(header_line, 2, 2), 'PG') THEN 'PG'
                    ELSE 'CO'
                END as header_type
            FROM read_csv_auto(
                '{pipe_path}', 
                header=FALSE, 
                columns={{'header_line':'VARCHAR'}},
                -- Increase buffer size for better throughput
                buffer_size=4194304,
                ignore_errors=true
            )
            """

            # Execute the query to create and populate the header table
            con.execute(query)

            # Wait for thread to finish
            thread.join()

            # Count and return the number of header lines
            header_count = con.execute("SELECT COUNT(*) FROM header").fetchone()[0]
            log.info(
                f"Processed {header_count} header lines into header table using named pipe in {'/dev/shm' if '/dev/shm' in pipe_path else 'tmp'}"
            )

            return header_count
        else:
            # Fall back to StringIO method if pipes cannot be created
            header_lines_str = "\n".join(header_lines)
            import io

            f = io.StringIO(header_lines_str)
            con.register("__HEADER_LINES__", f)

            query = """
            CREATE TABLE header AS
            SELECT 
                header_line,
                CASE
                    WHEN SUBSTRING(header_line, 2, 2) = 'HD' THEN 'HD'
                    WHEN SUBSTRING(header_line, 2, 2) = 'SQ' THEN 'SQ'
                    WHEN SUBSTRING(header_line, 2, 2) = 'RG' THEN 'RG'
                    WHEN SUBSTRING(header_line, 2, 2) = 'PG' THEN 'PG'
                    ELSE 'CO'
                END as header_type
            FROM read_csv_auto('__HEADER_LINES__', header=FALSE, columns={'header_line':'VARCHAR'})
            """

            con.execute(query)
            header_count = con.execute("SELECT COUNT(*) FROM header").fetchone()[0]
            log.info(
                f"Processed {header_count} header lines into header table using StringIO (fallback method)"
            )
            con.unregister("__HEADER_LINES__")
            return header_count

    except Exception as e:
        log.error(f"Error creating header tables: {e}")
        return 0
    finally:
        # Clean up the named pipe if it exists
        if use_pipe and pipe_path and os.path.exists(pipe_path):
            try:
                os.unlink(pipe_path)
                log.debug(f"Removed header named pipe: {pipe_path}")
            except Exception as e:
                log.warning(f"Could not remove header named pipe {pipe_path}: {e}")


def load_sam_file_optimized(
    con,
    file_path,
    match_reward=1,
    mismatch_penalty=-2,
    gap_open_penalty=5,
    gap_extension_penalty=2,
    keep_unused_references=False,
    threads=4,
    large_file_mode=True,
    use_shm=True,
):
    """Optimized SAM/BAM loading using direct pipe from samtools for BAM files or direct read for SAM."""
    original_db_threads = threads
    samtools_proc = None
    named_pipes = []  # Track all named pipes created

    try:
        # Configure DuckDB threads and optimizations
        current_max_threads = int(
            con.execute("SELECT current_setting('threads');").fetchone()[0]
        )
        log.info(f"DuckDB is configured with max {current_max_threads} threads.")
        effective_total_threads = min(original_db_threads, current_max_threads)

        # Performance optimizations
        con.execute("SET enable_progress_bar=true")
        con.execute("SET preserve_insertion_order=false")
        con.execute("SET default_null_order='nulls_first'")
        con.execute("SET force_compression='none'")
        con.execute("SET immediate_transaction_mode=true")

        # Determine file type
        base_ext, compression_type = get_file_extension(file_path)
        is_bam = base_ext == ".bam"
        is_sam = base_ext == ".sam"

        # Extract headers using samtools for any file type (BAM or SAM)
        header_lines = extract_sam_headers(file_path, threads=max(1, threads // 2))

        # Create header tables
        header_count = create_header_tables_from_raw_headers(con, header_lines)
        log.info(f"Created header table with {header_count} lines")

        # --- Only keep main/fast query paths below, remove fallbacks and extra checks ---

        if is_bam:
            log.info(f"Processing BAM file using samtools pipe: {file_path}")
            import os
            import uuid
            pipe_path = (
                f"/dev/shm/bam_filter_{uuid.uuid4()}"
                if os.path.exists("/dev/shm") and os.access("/dev/shm", os.W_OK)
                else os.path.join(tempfile.gettempdir(), f"bam_filter_{uuid.uuid4()}")
            )
            os.mkfifo(pipe_path)
            named_pipes.append(pipe_path)
            import threading

            def pipe_samtools_output():
                with open(pipe_path, "w") as pipe:
                    proc = subprocess.Popen(
                        ["samtools", "view"] + (["-@", str(max(1, threads // 2))] if threads > 1 else []) + [file_path],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        universal_newlines=True,
                        bufsize=-1,
                    )
                    for line in proc.stdout:
                        pipe.write(line)
                    proc.communicate()

            pipe_thread = threading.Thread(target=pipe_samtools_output)
            pipe_thread.daemon = True
            pipe_thread.start()

            load_query = f"""
            CREATE TEMP TABLE raw_input AS
            WITH parsed_input_with_raw_tags AS (
                SELECT
                    line,
                    fields,
                    CASE
                        WHEN array_length(fields) > 11
                        THEN array_to_string(list_slice(fields, 12, array_length(fields)), CHR(9))
                        ELSE NULL
                    END as _raw_tags_str
                FROM (
                    SELECT
                        line,
                        string_split(line, CHR(9)) AS fields
                    FROM read_csv(
                        '{pipe_path}',
                        delim='\\n', 
                        header=FALSE,
                        parallel=TRUE,
                        columns={{'line': 'VARCHAR'}},
                        buffer_size=268435456,
                        ignore_errors=true
                    )
                    WHERE NOT SUBSTRING(line, 1, 1) = '@'
                )
            )
            SELECT 
                qname,
                flag,
                rname,
                pos,
                mapq,
                cigar,
                rnext,
                pnext,
                tlen,
                seq,
                qual,
                raw_tags,
                query_length,
                NULLIF(REGEXP_EXTRACT(raw_tags, 'NM:i:([0-9]+)', 1), '')::INTEGER as nm_tag,
                NULLIF(REGEXP_EXTRACT(raw_tags, 'XO:i:([0-9]+)', 1), '')::INTEGER as xo_tag,
                NULLIF(REGEXP_EXTRACT(raw_tags, 'XG:i:([0-9]+)', 1), '')::INTEGER as xg_tag
            FROM parsed_input_with_raw_tags;
            """
            con.execute(load_query)
            log.info(f"Successfully loaded data using named pipe from samtools")
            pipe_thread.join()

        else:
            log.info(
                f"Processing {'compressed' if compression_type else 'uncompressed'} SAM file directly: {file_path}"
            )
            direct_load_query = f"""
            CREATE TEMP TABLE raw_input AS
            WITH raw_lines AS (
                SELECT
                    line
                FROM read_csv(
                    '{file_path}',
                    delim='\\n', 
                    header=FALSE,
                    parallel=TRUE,
                    comment='@',
                    columns={{'line': 'VARCHAR'}},
                    ignore_errors=true
                )
            ),
            extracted_fields AS (
                SELECT
                    line,
                    fields,
                    list_extract(fields, 1) as qname,
                    TRY_CAST(list_extract(fields, 2) AS SMALLINT) as flag,
                    list_extract(fields, 3) as rname,
                    TRY_CAST(list_extract(fields, 4) AS INTEGER) as pos,
                    TRY_CAST(list_extract(fields, 5) AS SMALLINT) as mapq,
                    NULLIF(list_extract(fields, 6), '*') as cigar,
                    list_extract(fields, 7) as rnext,
                    TRY_CAST(list_extract(fields, 8) AS INTEGER) as pnext,
                    TRY_CAST(list_extract(fields, 9) AS INTEGER) as tlen,
                    list_extract(fields, 10) as seq,
                    NULLIF(list_extract(fields, 11), '*') as qual,
                    CASE
                        WHEN array_length(fields) > 11
                        THEN array_to_string(list_slice(fields, 12, array_length(fields)), CHR(9))
                        ELSE NULL
                    END as raw_tags,
                    LENGTH(COALESCE(list_extract(fields, 10), '')) as query_length
                FROM (
                    SELECT
                        line,
                        string_split(line, CHR(9)) AS fields
                    FROM raw_lines
                )
            )
            SELECT 
                qname,
                flag,
                rname,
                pos,
                mapq,
                cigar,
                rnext,
                pnext,
                tlen,
                seq,
                qual,
                raw_tags,
                query_length,
                REGEXP_EXTRACT(raw_tags, 'NM:i:([0-9]+)', 1)::INTEGER as nm_tag,
                REGEXP_EXTRACT(raw_tags, 'XO:i:([0-9]+)', 1)::INTEGER as xo_tag,
                REGEXP_EXTRACT(raw_tags, 'XG:i:([0-9]+)', 1)::INTEGER as xg_tag
            FROM extracted_fields;
            """
            con.execute(direct_load_query)
            log.info("Successfully loaded SAM data using direct DuckDB read_csv")

        # --- Pre-calculate reference alignment counts (expensive operation done once) ---
        log.info("Calculating reference alignment counts...")
        con.execute(
            """
            CREATE TEMP TABLE reference_alignment_counts AS
            SELECT
                rname,
                COUNT(*) as aln_count
            FROM raw_input
            GROUP BY rname;
            """
        )

        # Optimize the reads query with a hash-based GROUP BY approach
        log.info("Creating reads table with optimized hash-based approach...")
        
        # Streamlined approach: Create reads table directly with window function
        # No need for temporary table or validation queries
        reads_query = """
        CREATE TABLE reads AS
        SELECT
            qname,
            flag,
            seq,
            qual,
            query_length AS length,
            rowid AS read_id
        FROM raw_input
        WHERE rowid IN (
            SELECT MIN(rowid)
            FROM raw_input
            GROUP BY qname
        );
        """
        con.execute(reads_query)
        
        # Simple log without expensive count query
        log.info("Created reads table using optimized window function approach")

        # --- Computing partition assignments using size-based approach ---
        log.info("Computing size-based partition assignments (256MB blocks, up to 100 partitions)...")
        
        # Calculate estimated row size in bytes based on alignment table column types
        row_size_bytes_query = """
        WITH sample_data AS (
            SELECT * FROM raw_input LIMIT 100
        ),
        column_sizes AS (
            SELECT
                -- Fixed-size columns (using standard sizes)
                4 AS read_id_size,          -- INTEGER (4 bytes)
                4 AS ref_id_size,           -- INTEGER (4 bytes)
                2 AS flag_size,             -- SMALLINT (2 bytes)
                4 AS pos_size,              -- INTEGER (4 bytes)
                2 AS mapq_size,             -- SMALLINT (2 bytes)
                4 AS cigar_id_size,         -- INTEGER (4 bytes)
                4 AS pnext_size,            -- INTEGER (4 bytes)
                4 AS tlen_size,             -- INTEGER (4 bytes)
                4 AS partition_id_size,     -- INTEGER (4 bytes)
                4 AS tag_ZA_size,           -- FLOAT (4 bytes)
                4 AS tag_ZS_size,           -- FLOAT (4 bytes)
                4 AS tag_id_size,           -- INTEGER (4 bytes)
                
                -- Variable-size column
                AVG(LENGTH(rnext)) AS avg_rnext_size
            FROM sample_data
        )
        SELECT 
            -- Sum all fixed sizes
            (read_id_size + ref_id_size + flag_size + pos_size + mapq_size + 
             cigar_id_size + pnext_size + tlen_size + partition_id_size + 
             tag_ZA_size + tag_ZS_size + tag_id_size +
             -- Add variable column with padding
             COALESCE(avg_rnext_size, 4) + 4) AS estimated_row_bytes
        FROM column_sizes;
        """
        
        estimated_row_bytes = con.execute(row_size_bytes_query).fetchone()[0]
        log.info(f"Estimated alignment row size: {estimated_row_bytes} bytes")
        
        # Calculate how many rows fit in a 256MB partition (268,435,456 bytes)
        partition_size_bytes = 268435456  # 256MB in bytes
        rows_per_partition = max(1, int(partition_size_bytes / estimated_row_bytes))
        
        # Get total alignment count to estimate total partitions needed
        total_alns = con.execute("SELECT SUM(aln_count) FROM reference_alignment_counts").fetchone()[0]
        estimated_partitions_needed = (total_alns + rows_per_partition - 1) // rows_per_partition
        
        # Cap at max_partitions if needed
        max_partitions = 100
        actual_partitions = min(estimated_partitions_needed, max_partitions)
        
        log.info(f"Data requires approximately {estimated_partitions_needed:,} partitions of 256MB each")
        
        # If we need more than max_partitions, adjust rows_per_partition
        if estimated_partitions_needed > max_partitions:
            rows_per_partition = (total_alns + max_partitions - 1) // max_partitions
            log.info(f"Limiting to {max_partitions} partitions with ~{rows_per_partition:,} rows per partition")
        else:
            log.info(f"Using {actual_partitions} partitions with ~{rows_per_partition:,} rows per partition")
        
        # Assign partitions based on accumulated row counts with INTEGER partition IDs
        con.execute(
            f"""
            CREATE TEMP TABLE ref_partition_map AS
            WITH ref_counts AS (
                SELECT
                    rac.rname,
                    rac.aln_count,
                    -- Order by alignment counts (largest first) for optimal packing
                    ROW_NUMBER() OVER (ORDER BY rac.aln_count DESC) AS ref_order
                FROM reference_alignment_counts rac
            ),
            -- Generate running counts for partition assignment
            running_counts AS (
                SELECT
                    rname,
                    aln_count,
                    -- Calculate cumulative rows to assign partitions
                    SUM(aln_count) OVER (ORDER BY ref_order) AS running_count,
                    -- Calculate partition number using INTEGER DIVISION to avoid fractional values
                    -- This ensures partition numbers start at 1 and are whole numbers
                    1 + (SUM(aln_count) OVER (ORDER BY ref_order) - 1) // {rows_per_partition} AS partition_number
                FROM ref_counts
            )
            SELECT
                rname, 
                -- Cast to INTEGER to ensure clean partition IDs
                -- Create sequential partition IDs based on the actual number of partitions we need
                -- If partition_number > actual_partitions, assign to partition 0
                CAST(
                    CASE 
                        WHEN partition_number > {actual_partitions} THEN 0
                        ELSE partition_number - 1  -- Zero-based partition IDs (0 to actual_partitions-1)
                    END 
                AS INTEGER) AS partition_id,
                aln_count
            FROM running_counts;
            """
        )
        
        # Log partition statistics
        log.info("Partition statistics:")
        partition_stats = con.execute("""
            SELECT 
                partition_id, 
                COUNT(*) AS num_refs, 
                SUM(aln_count) AS total_alns,
                ROUND(SUM(aln_count) * 100.0 / (SELECT SUM(aln_count) FROM ref_partition_map), 2) AS pct_of_total
            FROM ref_partition_map
            GROUP BY partition_id
            ORDER BY partition_id ASC
            LIMIT 10
        """).fetchall()
        
        for stat in partition_stats:
            log.info(f"Partition {stat[0]}: {stat[1]} refs, {stat[2]:,} alns ({stat[3]}% of total)")

        # --- Creating refs table with partition info ---
        log.info("Creating refs table using size-based partition assignments...")
        # Use a two-step approach that clearly avoids ambiguous rowid references
        refs_query = f"""
            -- First create sq_headers with extracted data from header
            CREATE TEMP TABLE sq_headers AS
            SELECT
                regexp_extract(header_line, 'SN:([^\\t]+)', 1) as SN,
                CAST(regexp_extract(header_line, 'LN:([0-9]+)', 1) AS INTEGER) as LN
            FROM header 
            WHERE header_type = 'SQ';
            
            -- Now extract distinct filtered refs with joins
            CREATE TEMP TABLE filtered_refs AS
            WITH direct_refs AS (
                -- Extract reference data directly from SQ headers with optimized extraction
                SELECT 
                    regexp_extract(header_line, 'SN:([^\\t]+)', 1) as rname,
                    CAST(regexp_extract(header_line, 'LN:([0-9]+)', 1) AS INTEGER) as ref_length
                FROM header 
                WHERE header_type = 'SQ'
            ),
            ref_counts AS (
                -- Pre-aggregate counts once for performance
                SELECT 
                    rname, 
                    COUNT(*) as aln_count
                FROM raw_input
                GROUP BY rname
            )
            SELECT 
                dr.rname,
                dr.ref_length,
                COALESCE(rc.aln_count, 0) as aln_count,
                COALESCE(rpm.partition_id, 0) as partition_id
            FROM direct_refs dr
            LEFT JOIN ref_counts rc ON dr.rname = rc.rname
            LEFT JOIN ref_partition_map rpm ON dr.rname = rpm.rname
            WHERE COALESCE(rc.aln_count, 0) > 0 OR {keep_unused_references};
            
            -- Then create the final refs table with its own rowid
            CREATE TABLE refs AS
            SELECT 
                rowid as ref_id,
                rname,
                ref_length,
                aln_count,
                partition_id
            FROM filtered_refs;
            
            -- Clean up temp tables
            DROP TABLE sq_headers;
            DROP TABLE filtered_refs;
        """
        con.execute(refs_query)

        # Optimize distinct operations for cigars using rowid
        log.info("Creating optimized cigars table using DuckDB's built-in rowid...")
        con.execute(
            """
            -- First create a temporary table for unique cigars
            CREATE TEMP TABLE unique_cigars AS
            SELECT 
                cigar
            FROM raw_input 
            GROUP BY cigar
            ORDER BY cigar;
            
            -- Then create the final cigars table with rowid
            CREATE TABLE cigars AS 
            SELECT 
                rowid as cigar_id,
                cigar
            FROM unique_cigars;
            
            -- Clean up the temp table
            DROP TABLE unique_cigars;
            """
        )
        
        # OPTIMIZED TAG PROCESSING: Use specialized fingerprinting technique
        # This approach creates tag signatures with fewer comparisons and less memory overhead
        log.info("Creating optimized tags table using bitmap fingerprinting...")
        con.execute(
            """
            -- Create tags with optimized bloom filter-like fingerprinting
            CREATE TABLE tags AS
            WITH tag_extraction AS (
                SELECT DISTINCT raw_tags
                FROM raw_input  
                WHERE raw_tags IS NOT NULL
            ),
            tag_classification AS (
                SELECT
                    raw_tags,
                    hash(raw_tags) AS tag_hash,
                    LENGTH(raw_tags) AS tag_length,
                    -- Create a bitmap signature from the tag content
                    (CASE WHEN raw_tags LIKE '%NM:i:%' THEN 1 ELSE 0 END) |
                    (CASE WHEN raw_tags LIKE '%AS:i:%' THEN 2 ELSE 0 END) |
                    (CASE WHEN raw_tags LIKE '%MD:Z:%' THEN 8 ELSE 0 END) |
                    (CASE WHEN raw_tags LIKE '%XS:i:%' THEN 16 ELSE 0 END) |
                    (CASE WHEN raw_tags LIKE '%ZA:%' THEN 128 ELSE 0 END) |
                    (CASE WHEN raw_tags LIKE '%ZS:%' THEN 256 ELSE 0 END) AS tag_bitmap
                FROM tag_extraction
            )
            SELECT 
                ROW_NUMBER() OVER (ORDER BY tag_bitmap, tag_length, tag_hash) AS tag_id,
                raw_tags,
                tag_hash,
                tag_length,
                tag_bitmap
            FROM tag_classification
            -- Smart ordering for alignment with memory access patterns
            ORDER BY tag_bitmap, tag_length, tag_hash;

            """
        )
        
        # OPTIMIZED ALIGNMENT PROCESSING: High-performance direct batch processing
        log.info("Creating intermediate table with high-performance batch processing...")
        con.execute(
            """            
            -- Create a simple tag lookup mapping to avoid expensive hash joins
            CREATE TEMP TABLE tag_lookup AS
            SELECT 
                raw_tags, 
                tag_id
            FROM tags;
            
            -- Process raw input directly with optimized joins 
            CREATE TEMP TABLE raw_input_with_ids AS
            SELECT
                p.flag,
                p.pos,
                p.mapq,
                p.rnext,
                p.pnext,
                p.tlen,
                p.query_length,
                p.nm_tag,
                p.xo_tag,
                p.xg_tag,
                r.read_id,
                rf.ref_id,
                rf.partition_id AS alignment_partition_id,
                c.cigar_id,
                -- Direct lookup with fallback
                COALESCE(tl.tag_id, 1) AS alignment_tag_id,
                -- Skip density class calculation which was expensive and rarely used
                CASE 
                    WHEN (SELECT COUNT(*) FROM raw_input WHERE rname = p.rname) > 1000000 THEN 'high_density'
                    WHEN (SELECT COUNT(*) FROM raw_input WHERE rname = p.rname) > 100000 THEN 'medium_density'  
                    ELSE 'low_density'
                END AS density_class
            FROM raw_input p
            JOIN reads r ON p.qname = r.qname
            JOIN refs rf ON p.rname = rf.rname
            JOIN cigars c ON p.cigar = c.cigar
            LEFT JOIN tag_lookup tl ON p.raw_tags = tl.raw_tags
            ORDER BY rf.partition_id, rf.ref_id, p.pos;
            
            -- Drop temporary indexes to free memory
            DROP INDEX IF EXISTS idx_reads_qname;
            DROP INDEX IF EXISTS idx_refs_rname; 
            DROP INDEX IF EXISTS idx_cigars_cigar;
            DROP INDEX IF EXISTS idx_tag_lookup_raw_tags;
            DROP TABLE IF EXISTS tag_lookup;
            """
        )

        # OPTIMIZED ALIGNMENT EXPORT PREPARATION
        # Create the alignments table directly from pre-sorted data
        alignments_query = f"""
            -- Create alignments table from pre-sorted data
            CREATE TABLE alignments AS
            SELECT
                read_id,
                ref_id,
                flag,
                pos,
                mapq,
                cigar_id,
                rnext,
                pnext,
                tlen,
                alignment_partition_id AS partition_id,
                -- Calculate ANI
                CASE WHEN query_length > 0
                    THEN ((1.0 - (COALESCE(nm_tag,0)::FLOAT / query_length)) * 100)::FLOAT
                    ELSE NULL
                END AS tag_ZA,
                -- Calculate shifted score
                CASE WHEN query_length > 0 THEN
                    ROUND((
                        ((query_length - COALESCE(nm_tag,0)) * {match_reward}) +
                        (COALESCE(nm_tag,0) * {mismatch_penalty}) -
                        (COALESCE(xo_tag,0) * {gap_open_penalty}) -
                        (COALESCE(xg_tag,0) * {gap_extension_penalty}) + 1.0
                    ) / query_length::FLOAT, 4)
                ELSE NULL END AS tag_ZS,
                alignment_tag_id AS tag_id
            FROM raw_input_with_ids;
            
            -- Create statistics to help optimizer with future queries
            ANALYZE alignments;
        """
        con.execute(alignments_query)

        # Remove @SQ lines from header table as this info is in 'refs'
        log.info("Removing @SQ lines from header table to reduce redundancy in Parquet export.")
        con.execute("DELETE FROM header WHERE header_type = 'SQ';")
        
        # Clean up temporary tables
        con.execute("DROP TABLE IF EXISTS raw_input")
        con.execute("DROP TABLE IF EXISTS reference_alignment_counts")
        con.execute("DROP TABLE IF EXISTS raw_input_with_ids")
        con.execute("DROP TABLE IF EXISTS ref_partition_map")

        # Force garbage collection
        import gc

        gc.collect()

        return True

    except Exception as e:
        log.error(f"Failed during optimized data loading: {e}")
        raise

    finally:
        # Clean up the samtools process if it exists
        if samtools_proc is not None and hasattr(samtools_proc, "poll"):
            try:
                if samtools_proc.poll() is None:  # Process still running
                    samtools_proc.terminate()
                    try:
                        samtools_proc.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        samtools_proc.kill()
                log.debug("Cleaned up samtools subprocess")
            except Exception as proc_e:
                log.warning(f"Error cleaning up samtools process: {proc_e}")

        # Clean up any named pipes
        for pipe in named_pipes:
            try:
                if os.path.exists(pipe):
                    os.unlink(pipe)
                    log.debug(f"Removed named pipe: {pipe}")
            except Exception as e:
                log.warning(f"Could not remove named pipe {pipe}: {e}")

        # Reset progress bar setting and threads
        try:
            con.execute("SET enable_progress_bar=false")
        except Exception:
            pass
        try:
            log.debug(f"Ensuring DuckDB threads are reset to {original_db_threads}.")
            con.execute(f"SET threads = {original_db_threads}")
        except Exception as thread_reset_e:
            log.warning(f"Could not reset DuckDB threads: {thread_reset_e}")


# This is a wrapper function that maintains compatibility with existing code
def load_sam_file(
    con,
    file_path,
    match_reward=1,
    mismatch_penalty=-2,
    gap_open_penalty=5,
    gap_extension_penalty=2,
    keep_unused_references=False,
    threads=4,
    large_file_mode=True,  # Add parameter for large file optimization
    use_stdin=False,  # For compatibility with existing code
):
    """
    Wrapper around optimized loading function that maintains compatibility
    with existing code while using the new optimized path.
    """
    return load_sam_file_optimized(
        con,
        file_path,
        match_reward=match_reward,
        mismatch_penalty=mismatch_penalty,
        gap_open_penalty=gap_open_penalty,
        gap_extension_penalty=gap_extension_penalty,
        keep_unused_references=keep_unused_references,
        threads=threads,
        large_file_mode=large_file_mode,
        use_shm=True,  # Always try to use shared memory for better performance
    )
