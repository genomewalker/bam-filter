import os
import logging
import subprocess
import shutil
import duckdb
import psutil
from pathlib import Path
import gzip
import bz2
import sys  # Import sys
import io  # Add for stdin handling
from bam_filter.utils import (
    check_tmp_dir_exists,
    get_compression_type,
    is_valid_file,
)
from bam_filter.sam_utils_db import (
    load_sam_file,
    COMMON_TAGS,
    dna_revcomp,
    str_reverse,
    get_file_extension,  # Import get_file_extension from here
)
from bam_filter.db_manager import DatabaseManager
import tempfile
import json
import math
from bam_filter import __version__  # Import version
import uuid  # Import uuid for unique temp file names
import pyarrow.parquet as pq  # Add pyarrow import

log = logging.getLogger("my_logger")


# Helper dictionary for SAM tag type characters
# (Extend this if COMMON_TAGS uses other types like float 'f' or array 'B')
SAM_TAG_TYPE_CHARS = {
    "INTEGER": "i",
    "VARCHAR": "Z",
    "FLOAT": "f",  # Added FLOAT type for ANI and shifted score
    # "REAL": "f", # Example if you add float tags
}


def format_sam_optional_tags(*args):
    """
    Reconstructs the SAM optional tags string.
    Args are expected to be: common_tag1_val, common_tag2_val, ..., raw_non_common_tags_str
    The order of common_tag_vals must match COMMON_TAGS.items() iteration.
    """
    tag_parts = []

    # Common tags are all args except the last one
    common_tag_values = args[:-1]
    # The last arg is the raw_non_common_tags_str
    raw_non_common_tags_str = args[-1]

    idx = 0
    for tag_name, duckdb_type in COMMON_TAGS.items():
        val = common_tag_values[idx]
        if val is not None:
            sam_type_char = SAM_TAG_TYPE_CHARS.get(
                duckdb_type, "Z"
            )  # Default to Z for unknown
            tag_parts.append(f"{tag_name}:{sam_type_char}:{val}")
        idx += 1

    if raw_non_common_tags_str is not None and raw_non_common_tags_str != "":
        # raw_non_common_tags_str is already tab-separated if it contains multiple tags
        tag_parts.append(raw_non_common_tags_str)

    if not tag_parts:
        return None  # Return None if no tags, SQL will handle COALESCE to '' or similar
    return "\t".join(tag_parts)


def get_table_info(con, table_name):
    """Gets column names, types, and total rows for a table."""
    try:
        # 'refs' does not need special quoting like 'references' did.
        # However, quoting all table names passed to DESCRIBE is a safe practice.
        quoted_table_name = f'"{table_name}"'
        columns_info = con.execute(f"DESCRIBE {quoted_table_name};").fetchall()
        col_names = [col[0] for col in columns_info]
        col_types = [col[1] for col in columns_info]
        total_rows = con.execute(
            f"SELECT COUNT(*) FROM {quoted_table_name};"
        ).fetchone()[0]
        return col_names, col_types, total_rows
    except Exception as e:
        log.error(f"Failed to get info for table {table_name}: {e}")
        return [], [], 0


def check_samtools():
    """Check if samtools is available in PATH"""
    if not shutil.which("samtools"):
        log.error(
            "samtools not found in PATH. Please install samtools and make sure it's in your PATH."
        )


def pipe_samtools_to_duckdb(input_file, db_manager, args):
    """
    Use samtools view to pipe BAM/SAM data directly to DuckDB without temporary files.

    Args:
        input_file (str): Path to the input BAM/SAM file
        db_manager (DatabaseManager): DuckDB connection manager
        args: Command-line arguments

    Returns:
        bool: True if successful, False otherwise
    """
    # Check if samtools is available
    if not shutil.which("samtools"):
        log.error("samtools not found in PATH. Please install samtools.")
        return False

    try:
        input_path = Path(input_file)
        base_ext, compression_type = get_file_extension(input_file)

        # Build samtools command
        cmd = ["samtools", "view", "-h"]  # -h to include headers

        # Add additional samtools options
        if hasattr(args, "threads") and args.threads:
            cmd.extend(
                ["-@", str(max(1, args.threads // 2))]
            )  # Use half threads for samtools

        cmd.append(str(input_path.resolve()))

        cmd_str = " ".join(cmd)
        log.info(f"Using command: {cmd_str}")

        # Try to use the PROGRAM feature of DuckDB for direct process data reading
        try:
            # First check if PROGRAM feature is available by running a small test
            test_result = db_manager.con.execute(
                "SELECT COUNT(*) FROM read_csv_auto('echo test')"
            ).fetchone()
            # If we can read from echo, then PROGRAM feature should work
            log.info(
                "DuckDB PROGRAM feature available, using direct pipe from samtools."
            )

            # Construct the query using the PROGRAM feature
            # This query must be identical to the raw_input creation in sam_utils_db.load_sam_file
            # apart from the input source.
            raw_input_creation_query = (
                """
                CREATE TEMP TABLE raw_input AS
                WITH parsed_input_with_raw_tags AS (
                    SELECT
                        line,
                        is_header,
                        fields, -- Keep fields for first 11 columns
                        -- Extract raw tags string once
                        CASE
                            WHEN NOT is_header AND array_length(fields) > 11
                            THEN array_to_string(list_slice(fields, 12, array_length(fields)), '\t')
                            ELSE NULL
                        END as _raw_tags_str
                    FROM ( -- Sub-select to define 'fields' and 'is_header' first
                        SELECT
                            line,
                            string_split(line, '\t') AS fields,
                            SUBSTRING(line, 1, 1) = '@' as is_header
                        FROM read_csv(
                            '"""
                + cmd_str
                + """',
                            delim='\\0', 
                            header=FALSE,
                            parallel=TRUE,
                            columns={'line': 'VARCHAR'},
                            buffer_size=268435456,
                            ignore_errors=true,
                            PROGRAM=true
                        )
                    )
                )
                SELECT 
                    line,
                    is_header,
                    CASE WHEN NOT is_header THEN list_extract(fields, 1) ELSE NULL END as qname,
                    CASE WHEN NOT is_header THEN TRY_CAST(list_extract(fields, 2) AS SMALLINT) ELSE NULL END as flag,
                    CASE WHEN NOT is_header THEN list_extract(fields, 3) ELSE NULL END as rname,
                    CASE WHEN NOT is_header THEN TRY_CAST(list_extract(fields, 4) AS INTEGER) ELSE NULL END as pos,
                    CASE WHEN NOT is_header THEN TRY_CAST(list_extract(fields, 5) AS SMALLINT) ELSE NULL END as mapq,
                    CASE WHEN NOT is_header THEN NULLIF(list_extract(fields, 6), '*') ELSE NULL END as cigar,
                    CASE WHEN NOT is_header THEN list_extract(fields, 7) ELSE NULL END as rnext,
                    CASE WHEN NOT is_header THEN TRY_CAST(list_extract(fields, 8) AS INTEGER) ELSE NULL END as pnext,
                    CASE WHEN NOT is_header THEN TRY_CAST(list_extract(fields, 9) AS INTEGER) ELSE NULL END as tlen,
                    CASE WHEN NOT is_header THEN list_extract(fields, 10) ELSE NULL END as seq,
                    CASE WHEN NOT is_header THEN NULLIF(list_extract(fields, 11), '*') ELSE NULL END as qual,
                    _raw_tags_str as raw_tags,
                    CASE
                        WHEN NOT is_header THEN LENGTH(COALESCE(list_extract(fields, 10), ''))
                        ELSE NULL
                    END as query_length,
                    CASE
                        WHEN NOT is_header THEN NULLIF(REGEXP_EXTRACT(_raw_tags_str, 'NM:i:([0-9]+)', 1), '')::INTEGER
                        ELSE NULL
                    END as nm_tag,
                    CASE
                        WHEN NOT is_header THEN NULLIF(REGEXP_EXTRACT(_raw_tags_str, 'XO:i:([0-9]+)', 1), '')::INTEGER
                        ELSE NULL
                    END as xo_tag,
                    CASE
                        WHEN NOT is_header THEN NULLIF(REGEXP_EXTRACT(_raw_tags_str, 'XG:i:([0-9]+)', 1), '')::INTEGER
                        ELSE NULL
                    END as xg_tag
                FROM parsed_input_with_raw_tags;
            """
            )

            # Execute the query which will start samtools and read its output
            db_manager.con.execute(
                raw_input_creation_query
            )  # MODIFIED: Use the full query
            log.info(
                "Successfully created raw_input table via DuckDB PROGRAM feature from samtools output."
            )
            # FIXME: After creating raw_input, this function should also create
            # header, reads, refs, cigars, and alignments tables,
            # similar to what load_sam_file does, to be a fully functional replacement.
            # For now, it correctly creates raw_input.
            return True

        except Exception as program_error:
            log.warning(
                f"DuckDB PROGRAM feature not available or failed: {program_error}"
            )

            # Fallback to alternative process piping method
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=-1,  # Use system default buffering
            )

            # Create a temporary named pipe for communication
            tmp_dir = Path(db_manager.temp_dir)
            pipe_path = tmp_dir / f"samtools_pipe_{uuid.uuid4()}"

            try:
                # Unix-style named pipe creation
                os.mkfifo(pipe_path)
                log.info(f"Created named pipe at {pipe_path}")

                # Start a thread to read from samtools and write to the pipe
                import threading

                def pipe_data():
                    with open(pipe_path, "w") as pipe:
                        for line in process.stdout:
                            pipe.write(line)

                thread = threading.Thread(target=pipe_data)
                thread.daemon = True  # Set as daemon so it closes with the main process
                thread.start()

                # Have DuckDB read from the named pipe file
                log.info("Loading data from samtools through named pipe")
                load_sam_file(
                    db_manager.con,
                    str(pipe_path),
                    match_reward=args.match_reward,
                    mismatch_penalty=args.mismatch_penalty,
                    gap_open_penalty=args.gap_open_penalty,
                    gap_extension_penalty=args.gap_extension_penalty,
                    threads=args.threads,
                    keep_unused_references=args.keep_unused_references,
                    use_stdin=False,  # Using a pipe, not stdin
                )

                # Check process status and log any errors
                stderr_out, _ = process.communicate()
                if stderr_out:
                    log.debug(f"samtools stderr output: {stderr_out}")

                return_code = process.returncode
                if return_code != 0:
                    log.error(f"samtools process exited with code {return_code}")
                    return False

                log.info("Successfully loaded data through named pipe")
                return True

            finally:
                # Clean up the named pipe
                if os.path.exists(pipe_path):
                    os.unlink(pipe_path)
                    log.debug(f"Removed named pipe: {pipe_path}")

                # Ensure process is terminated
                if process.poll() is None:
                    process.terminate()
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        process.kill()

    except Exception as e:
        log.error(f"Error during samtools piping: {e}")
        return False


def convert_bam_to_parquet(input_file, output_path, args):
    """Convert BAM or SAM file (compressed or uncompressed) to Parquet format"""
    tmp_dir = check_tmp_dir_exists(args.tmp_dir)
    compression = args.compression
    input_path = Path(input_file)
    base_ext, compression_type = get_file_extension(input_file)

    log.info(f"Converting {input_file} to Parquet format")

    if base_ext not in {".bam", ".sam"}:
        log.error(
            f"Unsupported input file base type: {base_ext} (compression: {compression_type})"
        )
        exit(1)

    sam_file_path = str(input_path.resolve())
    log.info(f"Processing input file: {sam_file_path}")

    # Try to calculate file size to adjust memory
    try:
        file_size_gb = os.path.getsize(sam_file_path) / 1024**3
        # Dynamically increase memory if needed
        if file_size_gb > 10 and not args.memory_limit:
            log.info(
                f"Large input file detected ({file_size_gb:.2f} GB), adjusting memory settings"
            )
            # Estimate memory needed (rule of thumb: 2.5x file size minimum, slightly more conservative)
            needed_memory_gb = max(int(file_size_gb * 2.5), 8)
            # Add a safety cap relative to available memory if possible
            try:
                available_memory_gb = psutil.virtual_memory().available / (1024**3)
                # Cap at 80% of available memory to leave room for OS and other processes
                capped_needed_memory_gb = min(
                    needed_memory_gb, int(available_memory_gb * 0.8)
                )
                if capped_needed_memory_gb < needed_memory_gb:
                    log.warning(
                        f"Calculated needed memory ({needed_memory_gb}GB) exceeds 80% of available RAM ({available_memory_gb:.2f}GB). Capping at {capped_needed_memory_gb}GB."
                    )
                needed_memory_gb = capped_needed_memory_gb
            except Exception as mem_check_e:
                log.warning(
                    f"Could not check available system memory: {mem_check_e}. Using calculated value {needed_memory_gb}GB without capping."
                )

            memory_limit = f"{needed_memory_gb}GB"
            log.info(
                f"Dynamically setting memory limit to {memory_limit} for large file processing."
            )
        elif args.memory_limit:
            memory_limit = args.memory_limit
            log.info(f"Using user-provided memory limit: {memory_limit}")
        else:
            memory_limit = args.memory_limit  # Which is None or default
            log.info("Using default memory limit (or None).")

    except Exception as size_calc_e:
        log.warning(
            f"Could not determine file size or set dynamic memory: {size_calc_e}. Using provided or default memory settings."
        )
        memory_limit = args.memory_limit

    try:
        # Ensure a file-based database is used
        db_file_name = f"bam_filter_conversion_{uuid.uuid4()}.db"
        db_path = Path(tmp_dir.name) / db_file_name
        log.info(f"Using temporary file-based database: {db_path}")

        with DatabaseManager(
            database=str(db_path),  # Use file path instead of :memory:
            temp_dir=tmp_dir.name,
            threads=args.threads,
            enable_progress=True,
            memory_limit=memory_limit,
            max_memory_pct=args.max_memory_pct,
        ) as db:
            log.info("Loading SAM/BAM file into database...")

            # Register Python UDFs
            db.con.create_function("dna_revcomp", dna_revcomp, [str], str)
            db.con.create_function("str_reverse", str_reverse, [str], str)
            log.debug("Registered dna_revcomp and str_reverse UDFs.")

            # Determine whether to use direct pipe from samtools
            use_pipe = base_ext == ".bam" and shutil.which("samtools") is not None

            if use_pipe:
                log.info("Using samtools to pipe BAM data directly to DuckDB")
                pipe_success = pipe_samtools_to_duckdb(sam_file_path, db, args)
                if not pipe_success:
                    log.warning(
                        "Piping with samtools failed, falling back to standard method"
                    )
                    use_pipe = False

            # If not using pipe or pipe failed, use standard load method
            if not use_pipe:
                load_sam_file(
                    db.con,
                    sam_file_path,
                    match_reward=args.match_reward,
                    mismatch_penalty=args.mismatch_penalty,
                    gap_open_penalty=args.gap_open_penalty,
                    gap_extension_penalty=args.gap_extension_penalty,
                    threads=args.threads,
                    keep_unused_references=args.keep_unused_references,
                )

            # Force checkpoint after load to ensure all data is saved
            # Force garbage collection to free memory

            # Run cleanup to free temporary files
            db.cleanup_temp_files()

            header_count = db.execute("SELECT COUNT(*) FROM header").fetchone()[0]
            try:
                alignment_count_res = db.execute(
                    "SELECT COUNT(*) FROM alignments"
                ).fetchone()
                alignment_count = alignment_count_res[0] if alignment_count_res else 0
            except Exception as e:
                log.warning(
                    f"Could not get alignment count (table might be empty or error occurred): {e}"
                )
                alignment_count = 0

            log.info(
                f"Loaded {header_count} header lines and {alignment_count} alignments"
            )

            # Inspect table storage after loading
            if args.enable_profiling:  # Only if profiling is enabled
                db.inspect_all_managed_tables_storage()

            # Create output directory if it doesn't exist
            output_dir = output_path.parent / output_path.stem
            output_dir.mkdir(parents=True, exist_ok=True)

            # Define tables to export - include dictionary tables for complete data storage
            tables_to_export = [
                "alignments",
                "header",
                "reads",
                "refs",
                "cigars",
                "tags",
            ]
            compression_level = (
                int(args.compression_level)
                if hasattr(args, "compression_level")
                and args.compression_level is not None
                else 11  # Default to 11
            )

            # Export all tables using the optimized export function
            for table_name in tables_to_export:
                export_success = db.export_optimized_table(
                    table_name, str(output_dir), compression_level=compression_level
                )
                if not export_success:
                    log.error(f"Failed to export {table_name}, see logs for details")
                else:
                    log.info(
                        f"Successfully exported {table_name} to {output_dir}/{table_name}.parquet"
                    )

    finally:
        # Before cleanup, check if we should keep the database file
        if hasattr(args, 'keep_db') and args.keep_db:
            # Determine the destination path
            if hasattr(args, 'db_path') and args.db_path:
                final_db_path = Path(args.db_path)
            else:
                # Default path based on output directory
                output_dir = output_path.parent / output_path.stem
                final_db_path = output_dir / f"{output_path.stem}.db"
            
            # Create parent directory if needed
            final_db_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Check if destination file exists
            if final_db_path.exists():
                if hasattr(args, 'overwrite') and args.overwrite:
                    log.warning(f"Overwriting existing database file: {final_db_path}")
                    try:
                        final_db_path.unlink()
                    except Exception as e:
                        log.error(f"Failed to remove existing database file: {e}")
                        if isinstance(tmp_dir, tempfile.TemporaryDirectory):
                            tmp_dir.cleanup()
                        return
                else:
                    log.error(
                        f"Destination database file exists: {final_db_path}. Use --overwrite to replace it."
                    )
                    # Continue with cleanup without copying the DB
                    if isinstance(tmp_dir, tempfile.TemporaryDirectory):
                        tmp_dir.cleanup()
                    # Skip DB persistence but don't exit (let other operations complete)
                    return
            
            # Copy the database file if it exists
            if "db_path" in locals() and db_path.exists():
                try:
                    log.info(f"Saving DuckDB database to: {final_db_path}")
                    # Get file size for progress reporting
                    db_size_mb = os.path.getsize(db_path) / (1024*1024)
                    log.info(f"Database file size: {db_size_mb:.2f} MB")
                    
                    # Copy with progress indication for large files
                    if db_size_mb > 100:  # Only show progress for larger files
                        with open(db_path, 'rb') as src, open(final_db_path, 'wb') as dst:
                            copied = 0
                            total_size = os.path.getsize(db_path)
                            buffer_size = 1024 * 1024  # 1MB buffer
                            
                            while True:
                                buffer = src.read(buffer_size)
                                if not buffer:
                                    break
                                
                                dst.write(buffer)
                                copied += len(buffer)
                                
                                # Print progress every ~10%
                                if copied % (total_size // 10) < buffer_size:
                                    percent = (copied / total_size) * 100
                                    log.info(f"Database copy progress: {percent:.1f}% ({copied/(1024*1024):.2f} MB of {db_size_mb:.2f} MB)")
                        
                        log.info(f"Successfully saved DuckDB database to {final_db_path}")
                    else:
                        # For smaller files, use simple copy
                        shutil.copy2(db_path, final_db_path)
                        log.info(f"Successfully saved DuckDB database to {final_db_path}")
                    
                    log.info(f"You can query it directly with: duckdb {final_db_path}")
                    
                    # Add helpful query examples
                    log.info("Sample queries you can run:")
                    log.info("  - SELECT COUNT(*) FROM alignments;")
                    log.info("  - SELECT * FROM refs LIMIT 10;")
                    log.info("  - SELECT rname, COUNT(*) as count FROM alignments JOIN refs ON alignments.ref_id = refs.ref_id GROUP BY rname ORDER BY count DESC LIMIT 20;")
                    log.info("  - SELECT AVG(tag_ZA) as avg_ani FROM alignments WHERE tag_ZA IS NOT NULL;")
                except Exception as e:
                    log.error(f"Failed to save the database file: {e}")
            else:
                log.warning("Database file not found or couldn't be saved.")
        
        # Standard cleanup (if not keeping or after copying)
        if isinstance(tmp_dir, tempfile.TemporaryDirectory):
            tmp_dir.cleanup()
        # Only remove the temporary DB file if not keeping it or if it was successfully copied
        if "db_path" in locals() and db_path.exists() and (not hasattr(args, 'keep_db') or not args.keep_db or 
                                                        (hasattr(args, 'keep_db') and args.keep_db and 
                                                         "final_db_path" in locals() and final_db_path.exists())):
            try:
                db_path.unlink()
                log.debug(f"Removed temporary database file: {db_path}")
            except Exception as e:
                log.warning(f"Could not remove temporary database file {db_path}: {e}")
                
    log.info(f"Conversion complete. Output saved to {output_path.stem}*.parquet")


def convert_parquet_to_sam(input_path, output_path, args):
    """Convert Parquet file(s) back to SAM format"""
    tmp_dir = check_tmp_dir_exists(args.tmp_dir)
    # input_path is expected to be the directory containing the parquet files/dirs
    input_dir = input_path

    log.info(f"Converting Parquet files from {input_dir} to SAM")

    # Define the expected base names for tables (can be files or directories)
    required_table_bases = [
        input_dir / "alignments.parquet",
        input_dir / "header.parquet",
        input_dir / "reads.parquet",
        input_dir / "refs.parquet",
        input_dir / "cigars.parquet",  # Add cigars dictionary
        input_dir / "tags.parquet",  # Add tags dictionary
    ]

    # Check if each required base path exists as either a file or directory
    missing_paths = []
    existing_paths = []  # Store the actual paths found (file or directory)
    for base_path in required_table_bases:
        if base_path.exists():  # Checks if path exists (file or dir)
            existing_paths.append(base_path)
        else:
            missing_paths.append(str(base_path))
    if missing_paths:
        log.error(
            f"Missing required Parquet files or directories: {', '.join(missing_paths)}"
        )
        log.error(
            f"Expected base paths (files or directories) in directory: {input_dir}/"
        )
        exit(1)
    try:
        # Ensure a file-based database is used
        db_file_name = f"parquet_filter_conversion_{uuid.uuid4()}.db"
        db_path = Path(tmp_dir.name) / db_file_name
        log.info(f"Using temporary file-based database: {db_path}")
        with DatabaseManager(
            database=str(db_path),  # Use file path instead of :memory:
            temp_dir=tmp_dir.name,
            threads=args.threads,
            enable_progress=True,
            memory_limit=args.memory_limit,
        ) as db:
            log.info("Loading Parquet files/directories into database...")
            # Prepare arguments for the new UDF
            # The UDF will take all common tag columns + the raw_tags column
            udf_arg_types = [COMMON_TAGS[key] for key in COMMON_TAGS.keys()] + [
                "VARCHAR"
            ]  # VARCHAR for raw_tags
            db.con.create_function(
                "format_sam_optional_tags", format_sam_optional_tags, udf_arg_types, str
            )
            log.debug("Registered format_sam_optional_tags UDF.")
            # Load all tables from the found paths (files or directories)
            for table_path in existing_paths:
                # Determine table name from the path stem (e.g., 'alignments.parquet' -> 'alignments')
                table_name = table_path.stem.split(".", 1)[0]
                # 'refs' does not need special quoting as it's not a SQL keyword.
                quoted_name = table_name
                # Determine the correct path string for read_parquet
                # If it's a directory, use a glob pattern to read all .parquet files within it
                if table_path.is_dir():
                    parquet_source_path = f"{str(table_path)}/*.parquet"
                    log.debug(
                        f"Reading from directory using glob: {parquet_source_path}"
                    )
                else:  # It's a single file
                    parquet_source_path = str(table_path)
                    log.debug(f"Reading from single file: {parquet_source_path}")
                try:
                    # read_parquet works with single files, lists of files, or glob patterns
                    db.execute(
                        f"CREATE OR REPLACE TABLE {quoted_name} AS SELECT * FROM read_parquet('{parquet_source_path}')"
                    )
                    count = db.execute(
                        f"SELECT COUNT(*) FROM {quoted_name}"
                    ).fetchone()[0]
                    log.info(f"Loaded {count} rows from {table_path}")
                except Exception as e:
                    # Add more context to the error log
                    log.error(
                        f"Failed to load {table_path} (using source path '{parquet_source_path}'): {e}"
                    )
                    exit(1)
            # Create the SAM file
            log.info(f"Exporting data to SAM file: {output_path}")
            # Determine compression for COPY command
            output_ext = output_path.suffix
            compression_string = "none"
            opener = open
            open_mode = "w"  # Mode for header writing
            if output_ext.endswith(".gz"):
                compression_string = "gzip"
                opener = gzip.open
                open_mode = "wt"
            elif output_ext.endswith(".bz2"):
                compression_string = "bzip2"
                opener = bz2.open
                open_mode = "wt"
            try:
                # --- Adjust @HD line for SO:unsorted ---
                log.info("Adjusting @HD header line to ensure SO:unsorted...")
                hd_line_result = db.execute(
                    "SELECT header_line FROM header WHERE header_type = 'HD' LIMIT 1"
                ).fetchone()
                if hd_line_result:
                    hd_line = hd_line_result[0]
                    import re

                    # Check if SO tag exists
                    if re.search(r"\tSO:[^\t]+", hd_line):
                        # Replace existing SO tag
                        new_hd_line = re.sub(r"\tSO:[^\t]+", "\tSO:unsorted", hd_line)
                    else:
                        # Append SO:unsorted
                        new_hd_line = hd_line + "\tSO:unsorted"
                    if new_hd_line != hd_line:
                        db.execute(
                            "UPDATE header SET header_line = ? WHERE header_line = ?",
                            (new_hd_line, hd_line),
                        )
                        log.info(f"Updated @HD line to: {new_hd_line}")
                    else:
                        log.info("@HD line already contains SO:unsorted or equivalent.")
                else:
                    # No @HD line found, insert a new one
                    # Using VN:1.6 as a reasonable default
                    new_hd_line = "@HD\tVN:1.6\tSO:unsorted"
                    db.execute(
                        "INSERT INTO header (header_line, header_type) VALUES (?, ?)",
                        (new_hd_line, "HD"),
                    )
                    log.info(f"Inserted new @HD line: {new_hd_line}")
                # --- End @HD adjustment ---
                # --- Add @PG line for this conversion command ---
                log.info("Adding @PG line for the Parquet to SAM conversion command...")
                import datetime

                pg_id = (
                    f"filterBAM.convert.{__version__}.{uuid.uuid4()}"  # Make ID unique
                )
                pg_pn = "filterBAM"
                pg_vn = __version__
                pg_cl = " ".join(sys.argv)  # Capture current command line
                pg_time = datetime.datetime.now().isoformat()
                pg_line = (
                    f"@PG\tID:{pg_id}\tPN:{pg_pn}\tVN:{pg_vn}\tCL:{pg_cl}\tDT:{pg_time}"
                )
                # Create a temporary table for the new PG line instead of inserting
                db.execute("DROP TABLE IF EXISTS temp_new_pg_line;")
                db.execute(
                    """
                    CREATE TEMP TABLE temp_new_pg_line (header_line VARCHAR, header_type VARCHAR);
                    """
                )
                db.execute(
                    "INSERT INTO temp_new_pg_line (header_line, header_type) VALUES (?, ?)",
                    (pg_line, "PG"),
                )
                log.info(f"Created temporary table for new @PG line: {pg_line}")
                # --- End @PG addition ---
                # Step 1: Create a temporary table with sorted header lines
                log.info("Creating temporary table with sorted header lines...")
                db.execute("DROP TABLE IF EXISTS temp_header_lines;")
                db.execute(
                    """
                    CREATE TEMP TABLE temp_header_lines AS
                    SELECT header_line
                    FROM (
                        SELECT
                            header_line,
                            header_type,
                            CASE header_type
                                WHEN 'HD' THEN 0
                                WHEN 'SQ' THEN 1
                                WHEN 'RG' THEN 2
                                WHEN 'PG' THEN 3
                                WHEN 'CO' THEN 4
                                ELSE 5
                            END AS sort_key
                        FROM (
                            SELECT header_line, header_type FROM header
                            UNION ALL
                            SELECT header_line, header_type FROM temp_new_pg_line
                        ) AS combined_headers
                    ) AS headers_with_key
                    ORDER BY
                        headers_with_key.sort_key,
                        headers_with_key.header_line;
                    """
                )
                header_count = db.execute(
                    "SELECT COUNT(*) FROM temp_header_lines"
                ).fetchone()[0]
                log.info(f"Created temporary header table with {header_count} lines.")
                # Construct the list of common tag columns for the UDF call
                common_tag_columns_for_udf = ", ".join(
                    [f"a.tag_{key}" for key in COMMON_TAGS.keys()]
                )
                # Ensure we have indexes on the numeric ID columns for faster joins
                log.info("Creating indexes on read_id and ref_id for faster joins...")
                try:
                    db.execute(
                        "CREATE INDEX IF NOT EXISTS idx_reads_read_id ON reads(read_id)"
                    )
                    db.execute(
                        "CREATE INDEX IF NOT EXISTS idx_refs_ref_id ON refs(ref_id)"
                    )
                    db.execute(
                        "CREATE INDEX IF NOT EXISTS idx_alignments_read_id ON alignments(read_id)"
                    )
                    db.execute(
                        "CREATE INDEX IF NOT EXISTS idx_alignments_ref_id ON alignments(ref_id)"
                    )
                    log.info("Indexes created successfully")
                except Exception as index_e:
                    log.warning(
                        f"Could not create indexes: {index_e}. Joins may be slower."
                    )
                # Get the total alignment count before using it in log message
                total_alignments = db.execute(
                    "SELECT COUNT(*) FROM alignments"
                ).fetchone()[0]
                # Create query that combines headers and alignments with explicit ordering
                combined_sql = f"""
                SELECT sam_output_line FROM (
                    SELECT header_line AS sam_output_line, 
                        0 AS section,
                        ROW_NUMBER() OVER() AS header_order
                    FROM temp_header_lines
                    UNION ALL
                    SELECT 
                        rd.qname || CHR(9) ||
                        a.flag::VARCHAR || CHR(9) ||
                        COALESCE(rf.rname, '*') || CHR(9) ||
                        COALESCE(a.pos::VARCHAR, '0') || CHR(9) ||
                        COALESCE(a.mapq::VARCHAR, '0') || CHR(9) ||
                        COALESCE(c.cigar, '*') || CHR(9) || 
                        COALESCE(a.rnext, '*') || CHR(9) ||
                        COALESCE(a.pnext::VARCHAR, '0') || CHR(9) ||
                        COALESCE(a.tlen::VARCHAR, '0') || CHR(9) ||
                        COALESCE(
                            CASE WHEN (rd.flag & 16) > 0 THEN 
                                dna_revcomp(rd.seq)
                            ELSE rd.seq END, 
                            '*' 
                        ) || CHR(9) ||
                        COALESCE(
                            CASE WHEN (rd.flag & 16) > 0 THEN str_reverse(rd.qual) ELSE rd.qual END, 
                            '*'
                        ) || CHR(9) ||
                        COALESCE(
                            format_sam_optional_tags({common_tag_columns_for_udf}, t.raw_tags), 
                            ''
                        ) AS sam_output_line,
                        1 AS section,
                        NULL AS header_order
                    FROM 
                        alignments a
                        JOIN reads rd ON a.read_id = rd.read_id
                        JOIN refs rf ON a.ref_id = rf.ref_id
                        JOIN cigars c ON a.cigar_id = c.cigar_id
                        JOIN tags t ON a.tag_id = t.tag_id
                ) combined_output
                ORDER BY 
                    section,
                    header_order
                """
                # Execute a single COPY command
                copy_command = f"""
                COPY ({combined_sql}) 
                TO '{str(output_path)}' 
                (FORMAT CSV, HEADER FALSE, DELIMITER E'\\x02', QUOTE '', ESCAPE '', COMPRESSION {compression_string})
                """
                log.info(f"Executing single COPY command with header ordering only...")
                db.execute(copy_command)
                log.info(
                    f"Successfully wrote {header_count} header lines and {total_alignments} unsorted alignment lines using single-query approach."
                )
                # Step 4: Clean up the temporary header tables
                db.execute("DROP TABLE temp_header_lines;")
                db.execute("DROP TABLE temp_new_pg_line;")
                log.debug("Dropped temporary header lines tables.")
            except Exception as e:
                log.error(f"Error exporting to SAM using combined COPY command: {e}")
                # Ensure temp tables are dropped even on error
                try:
                    db.execute("DROP TABLE IF EXISTS temp_header_lines;")
                    db.execute("DROP TABLE IF EXISTS temp_new_pg_line;")
                except Exception as drop_e:
                    log.warning(
                        f"Could not drop temp tables during error handling: {drop_e}"
                    )
                raise  # Re-raise the original error
    finally:
        if isinstance(tmp_dir, tempfile.TemporaryDirectory):
            tmp_dir.cleanup()
        # Clean up the specific database file if it exists
        if "db_path" in locals() and db_path.exists():
            try:
                db_path.unlink()
                log.debug(f"Removed temporary database file: {db_path}")
            except Exception as e:
                log.warning(f"Could not remove temporary database file {db_path}: {e}")


def read_file_list(file_list_path):
    """Read a list of files from a text file, one file path per line."""
    try:
        with open(file_list_path, "r") as f:
            files = [
                line.strip() for line in f if line.strip() and not line.startswith("#")
            ]
        # Validate all files exist
        missing = [f for f in files if not Path(f).exists()]
        if missing:
            log.error(
                f"The following files from the list do not exist: {', '.join(missing[:5])}"
                + (f" and {len(missing)-5} more..." if len(missing) > 5 else "")
            )
            exit(1)
        return files
    except Exception as e:
        log.error(f"Error reading file list: {e}")
        exit(1)


def convert_files(args):
    """Convert between SAM/BAM and Parquet formats"""
    if args.file_list:
        file_list = read_file_list(args.file_list)
        if args.output:
            output_path = Path(args.output)
            if output_path.suffix.lower() != ".parquet":
                log.warning(
                    f"Output path '{output_path}' does not end with .parquet. Using it as base name for output directory."
                )
            else:
                log.info(f"Output directory base name set to: {output_path}")
        else:
            first_file = Path(file_list[0])
            output_path = first_file.with_name(f"{first_file.stem}_combined")
            log.info(f"Output directory base name set to: {output_path}")
        output_dir = output_path.parent / output_path.stem
        expected_outputs = [
            output_dir / "alignments.parquet",
            output_dir / "header.parquet",
            output_dir / "reads.parquet",
            output_dir / "refs.parquet",
        ]
        output_exists = output_dir.exists() or any(f.exists() for f in expected_outputs)
        if output_exists:
            if not args.overwrite:
                log.error(
                    f"Output directory '{output_dir}' or its contents already exist. Use --overwrite to replace."
                )
                exit(1)
            else:
                log.warning(
                    f"Output directory '{output_dir}' or its contents already exist and will be overwritten (--overwrite specified)."
                )
                if output_dir.exists():
                    shutil.rmtree(output_dir)
                output_dir.parent.mkdir(parents=True, exist_ok=True)
        convert_file_list(file_list, output_path, args)
        return
    # --- Handle single file conversion ---
    input_path = Path(args.input)
    base_ext = None
    compression_type = None
    if not input_path.exists():
        log.error(f"Input path not found: {args.input}")
        exit(1)
    if input_path.is_dir():
        # Check if it looks like our Parquet output directory
        # Check if *any* file starting with alignments/header/reads exists (to handle multi-part files)
        has_alignments = any(input_path.glob("alignments.parquet*"))
        has_header = any(input_path.glob("header.parquet*"))
        has_reads = any(input_path.glob("reads.parquet*"))
        has_cigars = any(
            input_path.glob("cigars.parquet*")
        )  # Check for cigars dictionary
        has_tags = any(input_path.glob("tags.parquet*"))  # Check for tags dictionary

        if has_alignments and has_header and has_reads and has_cigars and has_tags:
            is_parquet_input = True
            parquet_input_dir = input_path
            base_ext = ".parquet"  # Treat directory as parquet input type
            log.info(
                f"Detected Parquet input directory with all required tables: {parquet_input_dir}"
            )
        else:
            missing_tables = []
            if not has_alignments:
                missing_tables.append("alignments.parquet*")
            if not has_header:
                missing_tables.append("header.parquet*")
            if not has_reads:
                missing_tables.append("reads.parquet*")
            if not has_cigars:
                missing_tables.append("cigars.parquet*")
            if not has_tags:
                missing_tables.append("tags.parquet*")

            log.error(
                f"Input directory '{input_path}' does not contain all expected Parquet files. Missing: {', '.join(missing_tables)}"
            )
            exit(1)
    elif input_path.is_file():
        base_ext, compression_type = get_file_extension(args.input)
        if base_ext == ".parquet":
            # Assume it's the 'marker' file and the directory is adjacent
            potential_dir = input_path.parent / input_path.stem
            has_alignments = any(potential_dir.glob("alignments.parquet*"))
            has_header = any(potential_dir.glob("header.parquet*"))
            has_reads = any(potential_dir.glob("reads.parquet*"))
            if potential_dir.is_dir() and has_alignments and has_header and has_reads:
                is_parquet_input = True
                parquet_input_dir = potential_dir
                log.info(
                    f"Detected Parquet input via marker file '{input_path}', using directory: {parquet_input_dir}"
                )
            else:
                log.error(
                    f"Input file '{input_path}' seems to be Parquet, but the corresponding directory '{potential_dir}' with required files was not found or incomplete."
                )
                exit(1)
        elif base_ext in {".bam", ".sam"}:
            is_parquet_input = False
            supported_compressions = {None, "gz", "bz2"}
            if compression_type not in supported_compressions:
                log.error(
                    f"Unsupported compression for {base_ext}: {compression_type}. Supported compressions are: none, gz, bz2"
                )
                exit(1)
        else:
            log.error(
                f"Unsupported input file format: {args.input}. Supported formats are SAM/BAM files or a Parquet directory."
            )
            exit(1)
    else:
        log.error(f"Input path is neither a file nor a directory: {args.input}")
        exit(1)

    # --- Determine output path and format ---
    if args.output:
        output_path = Path(args.output)
    else:
        # Auto-generate output path based on input type
        if is_parquet_input:
            # Input is Parquet dir, output should be SAM
            output_path = parquet_input_dir.with_suffix(".sam")
            log.info(f"Auto-generating output SAM file name: {output_path}")
        else:
            # Input is SAM/BAM, output should be Parquet dir marker
            if compression_type:
                # Handle cases like file.bam.gz -> file.parquet
                # Or file.sam.bz2 -> file.parquet
                stem = input_path.name.split(".", 1)[0]
                output_path = input_path.with_name(stem + ".parquet")
            else:
                output_path = input_path.with_suffix(".parquet")
            log.info(f"Auto-generating output Parquet marker file name: {output_path}")

    # --- Validate output format based on input format ---
    output_base_ext, output_compression = get_file_extension(str(output_path))
    if not is_parquet_input:  # Input is SAM/BAM
        if output_base_ext != ".parquet":
            log.error(
                f"Output format must be Parquet (directory) when input is BAM/SAM. Provided: {output_path}"
            )
            suggested_dir = output_path.with_suffix("")
            log.info(f"Output will be a directory named: {suggested_dir}/")
            log.info(
                "Files will be saved as alignments.parquet, header.parquet, reads.parquet, etc."
            )
            output_path = output_path.with_suffix(".parquet")
            log.info(f"Using marker file name: {output_path}")
        output_dir = output_path.parent / output_path.stem
        if output_dir.exists():
            if not args.overwrite:
                log.error(
                    f"Output directory '{output_dir}' already exists. Use --overwrite to replace its contents."
                )
                exit(1)
            else:
                log.warning(
                    f"Output directory '{output_dir}' already exists. Contents may be overwritten (--overwrite specified)."
                )
                # Deletion happens inside convert_bam_to_parquet if needed
        format_name = base_ext.upper()[1:]
        if compression_type:
            format_name = f"{format_name} ({compression_type})"
        log.info(f"Converting {format_name} to Parquet: {input_path} -> {output_dir}/")
        convert_bam_to_parquet(args.input, output_path, args)

    else:  # Input is Parquet directory
        if output_base_ext not in {
            ".sam",  # Removed BAM check
        }:
            log.warning(
                f"Output format should be SAM (e.g., ending in .sam, .sam.gz, .sam.bz2) when input is Parquet. Provided: {output_path}"
            )
            # Allow proceeding but warn user.
        # Check for SAM output specific compression support
        supported_output_compressions = {None, "gz", "bz2"}
        if output_compression not in supported_output_compressions:
            log.error(
                f"Unsupported output compression: {output_compression}. Use none, .gz, or .bz2 for SAM output."
            )
            exit(1)
        if output_path.exists():
            if not args.overwrite:
                log.error(
                    f"Output file '{output_path}' already exists. Use --overwrite to replace it."
                )
                exit(1)
            else:
                log.warning(
                    f"Output file '{output_path}' already exists and will be overwritten (--overwrite specified)."
                )
                try:
                    output_path.unlink()
                except OSError as e:
                    log.error(f"Could not remove existing output file: {e}")
                    exit(1)
        log.info(f"Converting Parquet to SAM: {parquet_input_dir} -> {output_path}")
        convert_parquet_to_sam(parquet_input_dir, output_path, args)


def convert_file_list(file_list, output_path, args):
    """Convert a list of files to the specified output format."""
    pass
