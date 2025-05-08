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
            db.con.execute("CHECKPOINT")

            # Force garbage collection to free memory
            import gc

            gc.collect()

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

            # Define tables to export
            tables_to_export = ["alignments", "header", "reads", "refs"]
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
        if isinstance(tmp_dir, tempfile.TemporaryDirectory):
            tmp_dir.cleanup()
        # Clean up the specific database file if it exists
        if "db_path" in locals() and db_path.exists():
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
    ]

    # Check if each required base path exists as either a file or directory
    missing_paths = []
    existing_paths = []  # Store the actual paths found (file or dir)
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

            # Register Python UDFs (needed if they are used in reconstruction logic, though less likely here)
            # It's good practice to register them if there's any chance SQL might call them.
            # The current SAM export SQL uses them.
            db.con.create_function("dna_revcomp", dna_revcomp, [str], str)
            db.con.create_function("str_reverse", str_reverse, [str], str)

            # Prepare arguments for the new UDF
            # The UDF will take all common tag columns + the raw_tags column
            udf_arg_types = [COMMON_TAGS[key] for key in COMMON_TAGS.keys()] + [
                "VARCHAR"
            ]  # VARCHAR for raw_tags
            db.con.create_function(
                "format_sam_optional_tags", format_sam_optional_tags, udf_arg_types, str
            )
            log.debug(
                "Registered dna_revcomp, str_reverse, and format_sam_optional_tags UDFs."
            )

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
                    # Keep \t here - Python regex on SAM line string
                    if re.search(r"\tSO:[^\t]+", hd_line):
                        # Replace existing SO tag
                        # Keep \t here - Python regex on SAM line string
                        new_hd_line = re.sub(r"\tSO:[^\t]+", "\tSO:unsorted", hd_line)
                    else:
                        # Append SO:unsorted
                        # Keep \t here - Python string concatenation for SAM line
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
                    # Keep \t here - Python string construction for SAM line
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
                # Keep \t here - Python string construction for SAM line
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
                    -- Explicitly create sort key column before ordering
                    SELECT header_line -- Only select the final column needed
                    FROM (
                        SELECT
                            header_line,
                            header_type,
                            CASE header_type -- Calculate sort key explicitly
                                WHEN 'HD' THEN 0
                                WHEN 'SQ' THEN 1
                                WHEN 'RG' THEN 2
                                WHEN 'PG' THEN 3
                                WHEN 'CO' THEN 4
                                ELSE 5
                            END AS sort_key
                        FROM (
                            -- Innermost UNION ALL combines original header and new PG line
                            SELECT header_line, header_type FROM header
                            UNION ALL
                            SELECT header_line, header_type FROM temp_new_pg_line
                        ) AS combined_headers
                    ) AS headers_with_key
                    ORDER BY
                        headers_with_key.sort_key, -- Order by the calculated key
                        headers_with_key.header_line; -- Secondary sort for stability
                    """
                )
                header_count = db.execute(
                    "SELECT COUNT(*) FROM temp_header_lines"
                ).fetchone()[0]
                log.info(f"Created temporary header table with {header_count} lines.")

                # NEW APPROACH: Use a single COPY command with ordering to ensure headers come first
                log.info(
                    "Using single-query approach with explicit ordering for SAM file creation..."
                )

                # Get the total alignment count before using it in log message
                total_alignments = db.execute(
                    "SELECT COUNT(*) FROM alignments"
                ).fetchone()[0]

                # Construct the list of common tag columns for the UDF call
                common_tag_columns_for_udf = ", ".join(
                    [f"a.tag_{key}" for key in COMMON_TAGS.keys()]
                )

                # Create query that combines headers and alignments with explicit ordering
                # Alignments will be output in an arbitrary (but likely insertion) order
                combined_sql = f"""
                SELECT sam_output_line FROM (
                    -- Headers (section 0) come first, include row_number for preserving order
                    SELECT 
                        header_line AS sam_output_line, 
                        0 AS section,
                        ROW_NUMBER() OVER() AS header_order -- Add row number to preserve original ordering
                    FROM temp_header_lines
                    
                    UNION ALL
                    
                    -- Alignments (section 1) come second
                    SELECT 
                        a.qname || CHR(9) || 
                        a.flag::VARCHAR || CHR(9) ||
                        COALESCE(a.rname, '*') || CHR(9) ||
                        COALESCE(a.pos::VARCHAR, '0') || CHR(9) ||
                        COALESCE(a.mapq::VARCHAR, '0') || CHR(9) ||
                        COALESCE(a.cigar, '*') || CHR(9) ||
                        COALESCE(a.rnext, '*') || CHR(9) ||
                        COALESCE(a.pnext::VARCHAR, '0') || CHR(9) ||
                        COALESCE(a.tlen::VARCHAR, '0') || CHR(9) ||
                        COALESCE(r.seq, '*') || CHR(9) ||
                        COALESCE(r.qual, '*') ||
                        COALESCE(
                            CHR(9) || format_sam_optional_tags({common_tag_columns_for_udf}, a.raw_tags), 
                            ''
                        ) AS sam_output_line,
                        1 AS section,
                        NULL AS header_order -- NULL for alignments as they're not sorted by this key
                    FROM 
                        alignments a
                        LEFT JOIN reads r ON a.qname = r.qname
                ) combined_output
                ORDER BY 
                    section, -- First by section (0=headers, 1=alignments)
                    header_order -- Then by original header order (preserves @HD, @SQ, @RG, @PG, @CO sequence)
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
                db.execute(
                    "DROP TABLE temp_new_pg_line;"
                )  # Drop the new temp table too
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

    # Determine if input is likely Parquet directory or SAM/BAM
    is_parquet_input = False
    parquet_input_dir = None
    base_ext = None
    compression_type = None

    if not input_path.exists():
        log.error(f"Input path not found: {args.input}")
        exit(1)

    if input_path.is_dir():
        # Check if it looks like our Parquet output directory
        expected_parquet_files = [
            input_path / "alignments.parquet",
            input_path / "header.parquet",
            input_path / "reads.parquet",
        ]
        # Check if *any* file starting with alignments/header/reads exists (to handle multi-part files)
        has_alignments = any(input_path.glob("alignments.parquet*"))
        has_header = any(input_path.glob("header.parquet*"))
        has_reads = any(input_path.glob("reads.parquet*"))

        if has_alignments and has_header and has_reads:
            is_parquet_input = True
            parquet_input_dir = input_path
            base_ext = ".parquet"  # Treat directory as parquet input type
            log.info(f"Detected Parquet input directory: {parquet_input_dir}")
        else:
            log.error(
                f"Input directory '{input_path}' does not contain the expected Parquet files (alignments.parquet*, header.parquet*, reads.parquet*)."
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
        # This case should ideally be caught by the initial exists() check
        log.error(f"Input path is neither a file nor a directory: {args.input}")
        exit(1)

    # --- Determine output path and format ---
    if args.output:
        output_path = Path(args.output)
    else:
        # Auto-generate output path based on input type
        if is_parquet_input:
            # Input is Parquet dir, output should be SAM
            # Use the directory name for the output SAM base name
            output_path = parquet_input_dir.with_suffix(".sam")
            log.info(f"Auto-generating output SAM file name: {output_path}")
        else:
            # Input is SAM/BAM, output should be Parquet dir marker
            if compression_type:
                # Handle cases like file.bam.gz -> file.parquet
                # Or file.bam -> file.parquet
                # Or file.sam.bz2 -> file.parquet
                # Find the first dot representing an extension
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
            # Adjust output_path to represent the directory marker file
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
            ".sam",
        }:  # Removed BAM check
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
    """Convert multiple BAM/SAM files from a list and combine them into one Parquet output set."""
    sample_file = file_list[0]
    base_ext, compression_type = get_file_extension(sample_file)
    if base_ext not in {".bam", ".sam"}:
        log.error("File list mode currently only supports BAM/SAM input files")
        exit(1)
    tmp_dir_obj = check_tmp_dir_exists(args.tmp_dir)  # Use the returned object
    tmp_dir_path = tmp_dir_obj.name  # Get the path string
    compression = args.compression
    log.info(f"Converting {len(file_list)} files to Parquet format")

    # Use a persistent temporary DB file for aggregation
    with tempfile.NamedTemporaryFile(
        delete=False, suffix=".db", dir=tmp_dir_path
    ) as tmp_db_file:
        agg_db_path = tmp_db_file.name
        log.info(f"Using aggregation database: {agg_db_path}")

    processed_files_count = 0
    temp_parquet_files = []  # Keep track of temp files for cleanup

    try:
        # Establish the main aggregation database connection
        with DatabaseManager(
            database=agg_db_path,
            temp_dir=tmp_dir_path,
            threads=args.threads,
            enable_progress=True,
            memory_limit=args.memory_limit,
            max_memory_pct=args.max_memory_pct,
        ) as db:
            # Create the schema in the aggregation database with normalized column names
            db.execute(
                """
                CREATE TABLE IF NOT EXISTS header (
                    header_line VARCHAR PRIMARY KEY,
                    header_type VARCHAR
                );
                
                CREATE TABLE IF NOT EXISTS refs (
                    rname VARCHAR PRIMARY KEY, -- Changed ref_name to rname
                    ref_length INTEGER,
                    aln_count INTEGER DEFAULT 0,
                    partition_id INTEGER DEFAULT 0
                );
                
                CREATE TABLE IF NOT EXISTS reads (
                    qname VARCHAR PRIMARY KEY,
                    seq VARCHAR,
                    qual VARCHAR
                );
                
                CREATE TABLE IF NOT EXISTS alignments (
                    qname VARCHAR,
                    flag SMALLINT,
                    rname VARCHAR,
                    pos INTEGER,
                    mapq SMALLINT,
                    cigar VARCHAR,
                    rnext VARCHAR,
                    pnext INTEGER,
                    tlen INTEGER,
                    raw_tags VARCHAR,
                    partition_id INTEGER DEFAULT 0
                );
            """
            )

            for i, file_path in enumerate(file_list):
                log.info(f"Processing file {i+1}/{len(file_list)}: {file_path}")

                # Generate unique names for temporary parquet files for this iteration
                run_uuid = uuid.uuid4()
                temp_header_pq = Path(tmp_dir_path) / f"temp_header_{run_uuid}.parquet"
                temp_refs_pq = Path(tmp_dir_path) / f"temp_refs_{run_uuid}.parquet"
                temp_reads_pq = Path(tmp_dir_path) / f"temp_reads_{run_uuid}.parquet"
                temp_align_pq = Path(tmp_dir_path) / f"temp_align_{run_uuid}.parquet"
                current_temp_files = [
                    temp_header_pq,
                    temp_refs_pq,
                    temp_reads_pq,
                    temp_align_pq,
                ]
                temp_parquet_files.extend(
                    current_temp_files
                )  # Add to overall list for final cleanup

                # Use a temporary file-based DB for loading each file
                with tempfile.NamedTemporaryFile(
                    delete=False, suffix=".db", dir=tmp_dir_path
                ) as tmp_load_db_file:
                    temp_load_db_path = tmp_load_db_file.name

                try:
                    with DatabaseManager(
                        database=temp_load_db_path,
                        temp_dir=tmp_dir_path,
                        threads=args.threads,
                        enable_progress=False,  # Less verbose logs
                        memory_limit=args.memory_limit,
                        max_memory_pct=args.max_memory_pct,
                    ) as temp_load_db:
                        log.info(
                            f"Loading file '{file_path}' into temporary database '{temp_load_db_path}'..."
                        )
                        load_status = load_sam_file(
                            temp_load_db.con,
                            file_path,
                            match_reward=args.match_reward,
                            mismatch_penalty=args.mismatch_penalty,
                            gap_open_penalty=args.gap_open_penalty,
                            gap_extension_penalty=args.gap_extension_penalty,
                            threads=args.threads,
                            keep_unused_references=args.keep_unused_references,
                        )
                        log.debug(f"Load status for {file_path}: {load_status}")

                        log.info(
                            f"Exporting data from '{file_path}' to temporary Parquet files..."
                        )
                        # Export tables to temporary Parquet files
                        temp_load_db.execute(
                            f"COPY header TO '{temp_header_pq}' (FORMAT PARQUET, COMPRESSION zstd, ROW_GROUP_SIZE 1000000, PER_THREAD_OUTPUT TRUE, PARQUET_VERSION v2, OVERWRITE TRUE)"
                        )
                        temp_load_db.execute(
                            f"COPY refs TO '{temp_refs_pq}' (FORMAT PARQUET, COMPRESSION zstd, ROW_GROUP_SIZE 1000000, PER_THREAD_OUTPUT TRUE, PARQUET_VERSION v2, OVERWRITE TRUE)"
                        )
                        temp_load_db.execute(
                            f"COPY reads TO '{temp_reads_pq}' (FORMAT PARQUET, COMPRESSION zstd, ROW_GROUP_SIZE 1000000, PER_THREAD_OUTPUT TRUE, PARQUET_VERSION v2, OVERWRITE TRUE)"
                        )
                        # Use optimized export for alignments if possible, otherwise standard COPY
                        optimized_export_succeeded = (
                            temp_load_db.export_optimized_alignments(
                                tmp_dir_path, compression_level=8
                            )
                        )

                        exported_align_path = (
                            Path(tmp_dir_path) / "alignments.parquet"
                        )  # Define expected path

                        if optimized_export_succeeded and exported_align_path.exists():
                            # If export_optimized_alignments succeeded and file exists, rename it
                            log.debug(
                                f"Renaming optimized export output: {exported_align_path} -> {temp_align_pq}"
                            )
                            try:
                                exported_align_path.rename(temp_align_pq)
                            except Exception as rename_e:
                                log.error(
                                    f"Failed to rename {exported_align_path} to {temp_align_pq}: {rename_e}"
                                )
                                log.warning(
                                    "Falling back to standard COPY TO Parquet for alignments."
                                )
                                # Attempt standard copy as fallback if rename fails
                                temp_load_db.execute(
                                    f"COPY alignments TO '{temp_align_pq}' (FORMAT PARQUET, COMPRESSION zstd, ROW_GROUP_SIZE 1000000, PER_THREAD_OUTPUT TRUE, PARQUET_VERSION v2, OVERWRITE TRUE)"
                                )
                        elif (
                            optimized_export_succeeded
                            and not exported_align_path.exists()
                        ):
                            # Handle case where export reported success but file is missing
                            log.error(
                                f"Optimized export reported success but expected file {exported_align_path} not found."
                            )
                            log.warning(
                                "Falling back to standard COPY TO Parquet for alignments."
                            )
                            temp_load_db.execute(
                                f"COPY alignments TO '{temp_align_pq}' (FORMAT PARQUET, COMPRESSION zstd, ROW_GROUP_SIZE 1000000, PER_THREAD_OUTPUT TRUE, PARQUET_VERSION v2, OVERWRITE TRUE)"
                            )
                        else:
                            # Fallback if optimized export fails initially
                            log.warning(
                                "Optimized alignment export failed or file missing, using standard COPY TO Parquet."
                            )
                            temp_load_db.execute(
                                f"COPY alignments TO '{temp_align_pq}' (FORMAT PARQUET, COMPRESSION zstd, ROW_GROUP_SIZE 1000000, PER_THREAD_OUTPUT TRUE, PARQUET_VERSION v2, OVERWRITE TRUE)"
                            )

                    # temp_load_db context manager handles closing and potential file deletion
                    log.info(
                        f"Merging data from '{file_path}' (via temp Parquet) into aggregation database '{agg_db_path}'..."
                    )

                    # Import from temporary Parquet files into the aggregation DB
                    db.execute(
                        f"INSERT OR IGNORE INTO header SELECT * FROM read_parquet('{temp_header_pq}')"
                    )
                    db.execute(
                        f"INSERT OR IGNORE INTO refs SELECT * FROM read_parquet('{temp_refs_pq}')"
                    )
                    db.execute(
                        f"INSERT OR IGNORE INTO reads SELECT * FROM read_parquet('{temp_reads_pq}')"
                    )
                    db.execute(
                        f"INSERT INTO alignments SELECT * FROM read_parquet('{temp_align_pq}')"
                    )

                    processed_files_count += 1
                    log.info(f"Finished merging data from '{file_path}'.")

                except Exception as e:
                    log.error(f"Error processing file {file_path}: {e}")
                    log.warning(f"Skipping file {file_path} due to error.")
                finally:
                    # Clean up the temporary database file for the individual load
                    if os.path.exists(temp_load_db_path):
                        try:
                            os.unlink(temp_load_db_path)
                        except Exception as e_unlink:
                            log.warning(
                                f"Could not remove temporary load database {temp_load_db_path}: {e_unlink}"
                            )
                    # Clean up the temporary parquet files for this iteration immediately
                    for pq_file in current_temp_files:
                        if pq_file.exists():
                            try:
                                os.unlink(pq_file)
                            except Exception as e_unlink_pq:
                                log.warning(
                                    f"Could not remove temporary parquet file {pq_file}: {e_unlink_pq}"
                                )

            if processed_files_count == 0:
                log.error("No files were processed successfully.")
                # Clean up aggregation DB file before exiting
                if os.path.exists(agg_db_path):
                    try:
                        os.unlink(agg_db_path)
                    except:
                        pass
                exit(1)

            # --- Cleanup unused references logic (remains the same) ---
            if not args.keep_unused_references:
                log.info(
                    "Cleaning up unused references and corresponding @SQ header lines from combined data..."
                )
            else:
                log.info(
                    "Keeping all references and header lines (--keep-unused-references specified)."
                )

            # --- Final counts and ANALYZE (remains the same) ---
            header_count = db.execute("SELECT COUNT(*) FROM header").fetchone()[0]
            alignment_count = db.execute("SELECT COUNT(*) FROM alignments").fetchone()[
                0
            ]
            read_count = db.execute("SELECT COUNT(*) FROM reads").fetchone()[0]
            ref_count = db.execute("SELECT COUNT(*) FROM refs").fetchone()[0]
            log.info(
                f"Combined database contains {header_count} header lines, {ref_count} references, {read_count} unique reads, and {alignment_count} alignments"
            )
            db.execute("ANALYZE alignments")
            db.execute("ANALYZE reads")
            db.execute("ANALYZE refs")
            db.execute("ANALYZE header")
            # --- End Final counts and ANALYZE ---

            # Inspect table storage after loading and merging all files
            if args.enable_profiling:  # Only if profiling is enabled
                db.inspect_all_managed_tables_storage()

            log.info("Writing combined Parquet output with optimized settings...")
            output_dir = output_path.parent / output_path.stem
            output_dir.mkdir(parents=True, exist_ok=True)

            # --- Export all tables using the optimized export function ---
            tables_to_export = ["alignments", "header", "reads", "refs"]
            compression_level = (
                int(args.compression_level)
                if hasattr(args, "compression_level")
                and args.compression_level is not None
                else 11  # Default to 11
            )

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
            # --- End Export logic ---

    finally:
        # Final cleanup of the aggregation database file
        if os.path.exists(agg_db_path):
            try:
                os.unlink(agg_db_path)
                log.debug(f"Removed aggregation database file: {agg_db_path}")
            except Exception as e:
                log.warning(
                    f"Could not remove aggregation database file {agg_db_path}: {e}"
                )

        # Cleanup any remaining temporary parquet files (should be empty if loop finished)
        for pq_file in temp_parquet_files:
            if pq_file.exists():
                try:
                    os.unlink(pq_file)
                except Exception as e_unlink_pq:
                    log.warning(
                        f"Could not remove leftover temporary parquet file {pq_file}: {e_unlink_pq}"
                    )

        # Cleanup the main temporary directory object if it was created by us
        if isinstance(tmp_dir_obj, tempfile.TemporaryDirectory):
            tmp_dir_obj.cleanup()

    log.info(f"Conversion complete. Combined output saved to {output_dir}")
