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
import signal  # Add signal import
from .parquet_utils import calculate_optimal_row_group_size, get_column_info, ColumnInfo
import pandas as pd  # Add at top with other imports
from .scoring import (
    calculate_ani_from_tags,
    calculate_alignment_score,
    calculate_shifted_scores,
)

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
    "MD": "VARCHAR",  # String encoding mismatched and deleted reference bases
    "PG": "VARCHAR",  # Program
    "XG": "INTEGER",  # Number of gap extensions
    "NM": "INTEGER",  # Edit distance
    "XM": "INTEGER",  # Number of mismatches
    "XN": "INTEGER",  # Number of ambiguous bases in reference
    "XO": "INTEGER",  # Number of gap opens
    "AS": "INTEGER",  # Alignment score
    "XS": "INTEGER",  # Suboptimal alignment score
    "YT": "VARCHAR",  # Tag indicating the type of read
    "RG": "VARCHAR",  # Read group
    "ZA": "FLOAT",  # ANI value (Average Nucleotide Identity) - Custom tag
    "ZS": "FLOAT",  # Shifted score - Custom tag
}

SAM_TAG_SCAN_SAMPLE_SIZE = (
    1_000_000  # Number of alignment lines to sample for tag scanning
)


def export_optimized_table(con, table_name, output_dir_str, compression_level=11):
    """
    Export a table to Parquet format with optimized settings.
    For the 'alignments' table, Hive partitioning by 'partition_id' is applied.

    Args:
        con: DuckDB connection object.
        table_name (str): Name of the table to export.
        output_dir_str (str): Path to the output base directory.
        compression_level (int): ZSTD compression level (1-22). (Note: This parameter is currently not used to form ZSTD_COMPRESSION_LEVEL or COMPRESSION_LEVEL in the SQL due to compatibility issues with some DuckDB versions. ZSTD default compression will be used.)

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
                export_path = (
                    output_dir / f"{table_name}.parquet"
                )  # Fallback to single file
        except Exception as e:
            log.warning(
                f"Could not describe table {table_name} to check for partition_id: {e}. Exporting without Hive partitioning to a single file."
            )
            export_path = (
                output_dir / f"{table_name}.parquet"
            )  # Fallback to single file
    else:
        export_path = output_dir / f"{table_name}.parquet"

    try:
        # Ensure the output directory for the table exists (especially for non-partitioned files)
        # For Hive partitioned, DuckDB creates the subdirectories.
        # For single files, their parent directory needs to exist.
        if (
            table_name != "alignments" or not partition_options
        ):  # If not alignments OR alignments but no partitioning
            export_path.parent.mkdir(parents=True, exist_ok=True)
        else:  # If alignments and partitioning is active, ensure the base output_dir exists
            output_dir.mkdir(parents=True, exist_ok=True)

        # Get column information for calculating row group size
        columns_info_list = []  # Default to empty list
        try:
            # Attempt single-argument call first, based on observed TypeError
            log.debug(
                f"Attempting get_column_info({table_name}) (single argument call)."
            )
            ret_val_single_arg = get_column_info(table_name)

            # Handle both a single ColumnInfo object and a list of ColumnInfo objects
            if isinstance(ret_val_single_arg, list):
                columns_info_list = ret_val_single_arg
            elif hasattr(
                ret_val_single_arg, "numpy_type"
            ):  # It's a single ColumnInfo object
                log.debug(
                    f"Converting single ColumnInfo object to list for '{table_name}'"
                )
                columns_info_list = [ret_val_single_arg]
            else:
                log.warning(
                    f"get_column_info(table_name) for '{table_name}' returned unexpected type {type(ret_val_single_arg)}. Value (first 100 chars): '{str(ret_val_single_arg)[:100]}'"
                )
                columns_info_list = []  # Ensure it's a list
        except TypeError as te_single_arg:
            log.warning(
                f"TypeError calling get_column_info(table_name) for '{table_name}': {te_single_arg}. This might indicate an outdated function signature or that it requires the connection. Trying get_column_info(con, table_name) as a fallback."
            )
            try:
                log.debug(
                    f"Attempting get_column_info(con, {table_name}) (two argument call as fallback)."
                )
                ret_val_two_args = get_column_info(con, table_name)
                if isinstance(ret_val_two_args, list):
                    columns_info_list = ret_val_two_args
                else:
                    log.warning(
                        f"get_column_info(con, table_name) for '{table_name}' returned type {type(ret_val_two_args)} (expected list). Value (first 100 chars): '{str(ret_val_two_args)[:100]}'"
                    )
                    columns_info_list = []  # Ensure it's a list for subsequent checks
            except Exception as e_two_args:
                log.warning(
                    f"Fallback call to get_column_info(con, table_name) also failed for '{table_name}': {e_two_args}"
                )
                columns_info_list = []  # Ensure it's a list
        except Exception as e_col_info:
            log.warning(
                f"An unexpected error occurred while calling get_column_info for '{table_name}': {e_col_info}"
            )
            columns_info_list = []  # Ensure it's a list

        row_group_size_option = ""
        if not columns_info_list:  # This covers empty list or if it was reset above
            log.warning(
                f"Could not retrieve valid column information for table '{table_name}'. Proceeding with default row group size."
            )
        else:
            try:
                # Ensure all elements in columns_info_list are not strings if calculate_optimal_row_group_size expects objects
                row_group_size = calculate_optimal_row_group_size(
                    con, table_name, columns_info_list
                )
                row_group_size_option = f", ROW_GROUP_SIZE {row_group_size}"
            except AttributeError as ae:
                log.error(
                    f"AttributeError during row group size calculation for '{table_name}': {ae}. This often means column_info_list has an unexpected structure. Skipping row group size optimization."
                )
            except Exception as e_calc_rg:
                log.error(
                    f"Error calculating row group size for '{table_name}': {e_calc_rg}. Skipping row group size optimization."
                )

        # Ensure no explicit compression level option (ZSTD_COMPRESSION_LEVEL or COMPRESSION_LEVEL) is used.
        # COMPRESSION ZSTD will use DuckDB's default ZSTD compression level.
        copy_query = f"""
        COPY (SELECT * FROM "{table_name}") 
        TO '{str(export_path)}' 
        (
            FORMAT PARQUET, 
            COMPRESSION ZSTD{row_group_size_option}{partition_options},
            OVERWRITE_OR_IGNORE TRUE
        );
        """
        log.info(
            f"Exporting table '{table_name}' to '{export_path}' with ZSTD compression (default level){row_group_size_option}{partition_options}..."
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


def scan_sam_tags(con):
    """Scan SAM file for all unique tags and log their distribution based on a sample."""
    log.info(
        f"Scanning for unique SAM tags (sampling up to {SAM_TAG_SCAN_SAMPLE_SIZE:,} alignments)..."
    )
    tags_query = f"""
    WITH sampled_raw_input AS (
        SELECT qname, _full_raw_tags  -- Changed from raw_tags to _full_raw_tags
        FROM raw_input
        WHERE NOT is_header AND _full_raw_tags IS NOT NULL -- Changed from raw_tags
        USING SAMPLE {SAM_TAG_SCAN_SAMPLE_SIZE} ROWS -- Reservoir sample
    ),
    tag_split AS (
        SELECT 
            regexp_split_to_table(sri._full_raw_tags, '\t') as tag, -- Changed from raw_tags
            sri.qname  -- Add qname to track per-alignment tags
        FROM sampled_raw_input sri
        -- _full_raw_tags IS NOT NULL is already handled by sampled_raw_input CTE
    ),
    tag_types AS (
        SELECT DISTINCT
            substring(tag, 1, 2) as tag_name,
            substring(tag, 4, 1) as tag_type,
            COUNT(*) as count, -- Number of occurrences of the tag in the sample
            COUNT(DISTINCT qname) as aln_count, -- Number of unique alignments in the sample having this tag
            -- Percentage of alignments in the *sample* that have this tag.
            (COUNT(DISTINCT qname) * 100.0 / NULLIF((SELECT COUNT(*) FROM sampled_raw_input), 0))::FLOAT as pct_alns
        FROM tag_split
        WHERE length(tag) >= 5
            AND position(':' in tag) > 0
        GROUP BY 1, 2
    )
    SELECT * FROM tag_types ORDER BY count DESC;
    """

    try:
        results = con.execute(tags_query).fetchall()
        log.info("\nTag Distribution:")
        if not results:
            log.info(
                "  No optional tags found in the input (all alignment records have 11 or fewer fields, or raw_tags field was empty/invalid)."
            )
            log.info(
                "  Individual tag columns (e.g., tag_NM, tag_AS) in the alignments table will be NULL."
            )
            return False  # Explicitly return False if no tags are found

        log.info(
            f"{'Tag':4} | {'Type':4} | {'Count':>10} | {'Alignments':>10} | {'Coverage %':>9} | {'In COMMON_TAGS':>13}"
        )
        log.info("-" * 70)
        for tag_name, tag_type, count, aln_count, pct_alns in results:
            is_common = tag_name in COMMON_TAGS
            log.info(
                f"{tag_name:4} | {tag_type:4} | {count:>10,d} | {aln_count:>10,d} | "
                f"{pct_alns:>8.2f}% | {'Yes' if is_common else 'No':^13}"
            )
        return True  # Return True if tags were found and processed
    except Exception as e:
        log.error(f"Error scanning tags: {e}")
        return False


def load_sam_file(
    con,
    file_path,
    match_reward=1,
    mismatch_penalty=-2,
    gap_open_penalty=5,
    gap_extension_penalty=2,
    keep_unused_references=False,
    threads=4,
):
    """Optimized SAM/BAM loading with direct schema definition and tab parsing."""
    original_db_threads = threads  # Store the total requested threads
    try:
        try:
            # Use SELECT current_setting for modern DuckDB
            current_max_threads = int(
                con.execute("SELECT current_setting('threads');").fetchone()[0]
            )
            log.info(f"DuckDB is configured with max {current_max_threads} threads.")
            # Use the lower of the requested threads and the DB's max threads as the effective total
            effective_total_threads = min(original_db_threads, current_max_threads)
            if effective_total_threads < original_db_threads:
                log.warning(
                    f"Requested {original_db_threads} threads, but DB max is {current_max_threads}. Using {effective_total_threads} total."
                )
            else:
                # Use the originally passed value if within limits or if DB max couldn't be determined
                effective_total_threads = original_db_threads
        except Exception:
            # Fallback for older DuckDB versions or if setting query fails
            log.warning(
                f"Could not query current thread setting. Using requested threads value ({original_db_threads}) directly for total."
            )
            effective_total_threads = original_db_threads

        # Ensure we have at least 1 thread total
        effective_total_threads = max(1, effective_total_threads)

        # Enable progress bar for this operation
        con.execute("SET enable_progress_bar=true")

        read_query = f"""
            -- Create temporary table to hold raw input with pre-split fields and raw_tags
            CREATE TEMP TABLE raw_input AS
            WITH split_lines AS (
                SELECT
                    column0 AS line,
                    string_split(column0, '\t') AS fields
                FROM read_csv_auto(
                    '{file_path}',
                    delim='\\0',
                    header=FALSE
                )
            )
            SELECT 
                line,
                fields, -- Store the full fields array
                SUBSTRING(line, 1, 1) = '@' as is_header,
                list_extract(fields, 1) as qname,
                TRY_CAST(list_extract(fields, 2) AS SMALLINT) as flag,
                list_extract(fields, 3) as rname,
                TRY_CAST(list_extract(fields, 4) AS INTEGER) as pos,
                TRY_CAST(list_extract(fields, 5) AS SMALLINT) as mapq,
                NULLIF(list_extract(fields, 6), '*') as cigar,
                NULLIF(list_extract(fields, 7), '*') as rnext,
                TRY_CAST(list_extract(fields, 8) AS INTEGER) as pnext,
                TRY_CAST(list_extract(fields, 9) AS INTEGER) as tlen,
                list_extract(fields, 10) as seq, 
                list_extract(fields, 11) as qual,
                CASE
                    WHEN array_length(fields) > 11 THEN
                        -- Keep tabs as the separator for the optional tags
                        array_to_string(list_slice(fields, 12, array_length(fields)), '\t')
                    ELSE NULL
                END as _full_raw_tags -- Construct from fields array
            FROM split_lines;

            -- Create header table first
            CREATE TABLE header AS
            SELECT 
                line as header_line,
                CASE
                    WHEN SUBSTRING(line, 2, 2) = 'HD' THEN 'HD'
                    WHEN SUBSTRING(line, 2, 2) = 'SQ' THEN 'SQ'
                    WHEN SUBSTRING(line, 2, 2) = 'RG' THEN 'RG'
                    WHEN SUBSTRING(line, 2, 2) = 'PG' THEN 'PG'
                    ELSE 'CO'
                END as header_type
            FROM raw_input
            WHERE is_header;
        """

        con.execute(read_query)

        # --- BEGIN ADDED DEBUG LOGGING ---
        if log.isEnabledFor(logging.DEBUG):
            debug_raw_input_query = """
            SELECT
                line, -- Keep line for full context
                qname,
                flag,
                rname,
                pos,
                mapq,
                cigar,
                _full_raw_tags, -- Check this value carefully (renamed from raw_tags)
                array_length(fields) as num_fields_from_array,
                list_slice(fields, 12, NULLIF(array_length(fields),0)) as optional_fields_from_array
            FROM raw_input
            WHERE NOT is_header AND line NOT LIKE '@%' AND length(line) > 10
            LIMIT 5;
            """
            log.debug("Sample of raw_input table content and field splitting:")
            try:
                sample_raw_input = con.execute(debug_raw_input_query).fetchall()
                if sample_raw_input:
                    for row_idx, row in enumerate(sample_raw_input):
                        log.debug(f"Raw_input Sample Row {row_idx + 1}:")
                        log.debug(f"  Line (first 100 chars): '{str(row[0])[:100]}...'")
                        log.debug(f"  QNAME (parsed): {row[1]}")
                        log.debug(f"  FLAG (parsed): {row[2]}")
                        log.debug(f"  RNAME (parsed): {row[3]}")
                        log.debug(
                            f"  Full Raw Tags (in raw_input): {row[7]}"
                        )  # Updated field name
                        log.debug(
                            f"  Num Fields (from stored 'fields' array): {row[8]}"
                        )
                        log.debug(
                            f"  Optional Fields (from stored 'fields' array, 12+): {row[9]}"
                        )
                else:
                    log.debug(
                        "No alignment-like lines found in raw_input for sampling, or raw_input is empty."
                    )
            except Exception as e_debug:
                log.debug(
                    f"Could not fetch sample from raw_input for debugging: {e_debug}"
                )
        # --- END ADDED DEBUG LOGGING ---

        # Now scan tags after raw_input and header tables are created
        detected_tags = scan_sam_tags(con)  # scan_sam_tags now uses _full_raw_tags
        # common_tags_extraction_sql will parse from _full_raw_tags by default
        common_tags_extraction_sql = create_tag_extract_sql()

        # Create remaining tables (refs, reads) - no change here
        tables_query = f"""
            -- Create refs table
            CREATE TABLE refs AS
            WITH sq_lines AS (
                SELECT 
                    header_line,
                    regexp_extract(header_line, 'SN:([^\\t]+)', 1) as rname, 
                    CAST(regexp_extract(header_line, 'LN:([0-9]+)', 1) AS INTEGER) as ref_length
                FROM header 
                WHERE header_type = 'SQ'
            ), 
            aln_counts AS (
                SELECT 
                    rname, -- Use pre-parsed rname from raw_input
                    COUNT(*) as cnt
                FROM raw_input 
                WHERE NOT is_header 
                    AND qname IS NOT NULL -- Use qname to filter valid alignment lines
                    AND rname != '*'     -- Use pre-parsed rname
                GROUP BY rname
            )
            SELECT 
                sq.rname, 
                sq.ref_length,
                COALESCE(ac.cnt, 0) as aln_count,
                0 as partition_id
            FROM sq_lines sq
            INNER JOIN aln_counts ac ON sq.rname = ac.rname 
            WHERE ac.cnt > 0;

            -- Create reads table
            CREATE TABLE reads AS
            WITH raw_reads AS (
                SELECT DISTINCT
                    qname,
                    seq,  -- Use pre-parsed seq from raw_input
                    qual, -- Use pre-parsed qual from raw_input
                    flag  -- Use pre-parsed flag from raw_input
                FROM raw_input 
                WHERE NOT is_header 
                    AND qname IS NOT NULL -- Use qname to filter valid alignment lines
            )
            SELECT DISTINCT 
                qname,
                CASE 
                    WHEN flag & 16 > 0 THEN dna_revcomp(CASE WHEN seq = '*' THEN NULL ELSE seq END)
                    ELSE CASE WHEN seq = '*' THEN NULL ELSE seq END
                END as seq,
                CASE 
                    WHEN flag & 16 > 0 THEN str_reverse(CASE WHEN qual = '*' THEN NULL ELSE qual END)
                    ELSE CASE WHEN qual = '*' THEN NULL ELSE qual END
                END as qual
            FROM raw_reads;
        """

        con.execute(tables_query)

        # Now create alignments table with partition stats
        # Use proper tab separator consistently throughout this code
        common_tag_keys_sql_list = ", ".join([f"'{key}'" for key in COMMON_TAGS.keys()])
        filtered_raw_tags_expression = f"""
            NULLIF(TRIM(BOTH '\t' FROM array_to_string(
                list_filter(
                    string_split(_full_raw_tags, '\t'),
                    tag_item -> array_length(string_split(tag_item, ':')) = 3 AND substring(tag_item, 1, 2) NOT IN ({common_tag_keys_sql_list})
                ),
                '\t'
            )), '') 
        """

        # Add alignment score and ANI calculation expressions with explicit casting
        # These calculations will be used internally but won't be included in the final table
        score_calc_expression = f"""
            -- Calculate alignment score components with explicit casting
            CAST(COALESCE(CAST(tag_NM AS INTEGER), 0) AS INTEGER) as _num_mismatches,
            CASE 
                WHEN cigar IS NOT NULL THEN LENGTH(REGEXP_REPLACE(cigar, '[0-9]', '')) - CAST(COALESCE(CAST(tag_NM AS INTEGER), 0) AS INTEGER)
                ELSE LENGTH(COALESCE(raw_input.seq, '')) - CAST(COALESCE(CAST(tag_NM AS INTEGER), 0) AS INTEGER)
            END as _num_matches,
            CAST(COALESCE(CAST(tag_XO AS INTEGER), 0) AS INTEGER) as _num_gaps,
            CAST(COALESCE(CAST(tag_XG AS INTEGER), 0) AS INTEGER) as _gap_extensions,
            
            -- Calculate alignment score with explicit casting
            (
                -- Match reward * matches
                (CASE WHEN cigar IS NOT NULL THEN LENGTH(REGEXP_REPLACE(cigar, '[0-9]', '')) - CAST(COALESCE(CAST(tag_NM AS INTEGER), 0) AS INTEGER)
                 ELSE LENGTH(COALESCE(raw_input.seq, '')) - CAST(COALESCE(CAST(tag_NM AS INTEGER), 0) AS INTEGER) END) * {match_reward}
            ) - (
                -- Mismatch penalty * mismatches
                CAST(COALESCE(CAST(tag_NM AS INTEGER), 0) AS INTEGER) * {mismatch_penalty}
            ) - (
                -- Gap open penalty * gap opens
                CAST(COALESCE(CAST(tag_XO AS INTEGER), 0) AS INTEGER) * {gap_open_penalty}
            ) - (
                -- Gap extension penalty * gap extensions
                CAST(COALESCE(CAST(tag_XG AS INTEGER), 0) AS INTEGER) * {gap_extension_penalty}
            ) as _alignment_score,
            
            -- Calculate ANI (Average Nucleotide Identity) with explicit casting
            CASE
                WHEN LENGTH(COALESCE(raw_input.seq, '')) > 0 THEN 
                    1.0 - (CAST(COALESCE(CAST(tag_NM AS INTEGER), 0) AS FLOAT) / LENGTH(COALESCE(raw_input.seq, '')))
                ELSE NULL
            END as _ani
        """

        alignments_query = f"""
            -- Create partition stats
            CREATE TEMP TABLE partition_stats AS
            WITH alignment_counts AS (
                SELECT 
                    rname, -- Use pre-parsed rname from raw_input
                    COUNT(*) as alignment_count
                FROM raw_input 
                WHERE NOT is_header 
                    AND qname IS NOT NULL -- Use qname to filter valid alignment lines
                    AND rname != '*'     -- Use pre-parsed rname
                GROUP BY rname
            ),
            ranked_refs AS (
                SELECT 
                    rname, 
                    alignment_count,
                    SUM(alignment_count) OVER () as total_alignments,
                    (alignment_count * 100.0 / SUM(alignment_count) OVER ())::FLOAT as pct_of_total,
                    SUM(alignment_count) OVER (ORDER BY alignment_count DESC) as cumulative_alns,
                    (SUM(alignment_count) OVER (ORDER BY alignment_count DESC) * 100.0 / 
                     SUM(alignment_count) OVER ())::FLOAT as cumulative_pct
                FROM alignment_counts
            )
            SELECT 
                rname as rname, 
                alignment_count as alignments,
                pct_of_total,
                total_alignments as total_alns,
                cumulative_pct,
                CASE
                    WHEN cumulative_pct <= 1.0 THEN 99
                    WHEN cumulative_pct <= 5.0 THEN 98
                    WHEN cumulative_pct <= 15.0 THEN 97
                    WHEN cumulative_pct <= 30.0 THEN 96
                    ELSE GREATEST(0, 95 - FLOOR((cumulative_pct - 30.0) * 95.0 / 70.0)::INTEGER)
                END as partition_id
            FROM ranked_refs;

            -- Create temporary calculation table with all the metrics
            CREATE TEMP TABLE temp_alignment_calculations AS 
            WITH parsed_alignments AS (
                SELECT 
                    raw_input.qname,    -- Explicitly specify raw_input.qname
                    raw_input.flag,     -- Directly from raw_input
                    raw_input.rname,    -- Directly from raw_input
                    raw_input.pos,      -- Directly from raw_input
                    raw_input.mapq,     -- Directly from raw_input
                    raw_input.cigar,    -- Directly from raw_input
                    raw_input.rnext,    -- Directly from raw_input
                    raw_input.pnext,    -- Directly from raw_input
                    raw_input.tlen,     -- Directly from raw_input
                    LENGTH(COALESCE(raw_input.seq, '')) as query_length, -- Add query length
                    raw_input._full_raw_tags, -- Keep full tags temporarily for parsing
                    {common_tags_extraction_sql}, -- Operates on _full_raw_tags with tab-based extraction
                    {filtered_raw_tags_expression} AS raw_tags, -- This is the new filtered raw_tags
                    {score_calc_expression} -- Add score calculations with fixed castings
                FROM raw_input 
                LEFT JOIN reads ON raw_input.qname = reads.qname -- Changed to LEFT JOIN to preserve all alignments
                WHERE NOT raw_input.is_header 
                    AND raw_input.qname IS NOT NULL -- Use qname to filter valid alignment lines
            )
            SELECT 
                qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, -- Standard fields
                query_length, -- Include query length
                raw_tags, -- The new filtered raw_tags
                {', '.join([f'tag_{key}' for key in COMMON_TAGS.keys()])}, -- All common tag columns
                _alignment_score, -- Keep temporarily for shifted score calculation
                _ani -- Keep temporarily for ANI tag
            FROM parsed_alignments;
            
            -- Now create the final alignments table without the intermediate calculation fields
            CREATE TABLE alignments AS
            SELECT 
                qname, flag, t.rname, pos, mapq, cigar, rnext, pnext, tlen, -- Standard fields
                query_length, -- Include query length
                raw_tags, -- The filtered raw_tags
                {', '.join([f'tag_{key}' for key in COMMON_TAGS.keys()])}, -- All common tag columns
                COALESCE(ps.partition_id, 0) as partition_id
            FROM temp_alignment_calculations t
            LEFT JOIN partition_stats ps ON t.rname = ps.rname;
        """

        con.execute(alignments_query)

        # After creating the alignments table, directly add ANI and shifted_score values as tags
        # Calculate from the temporary table that contains the needed metrics
        log.info("Calculating shifted scores and adding custom tags directly...")

        # Calculate shifted scores and add them as ZS and ZA tags
        shifted_scores_query = f"""
        -- First, get minimum score for normalization (run once for efficiency)
        WITH min_score AS (
            SELECT MIN(_alignment_score) as min_score
            FROM temp_alignment_calculations
        )
        
        -- Add the ZS (shifted_score) and ZA (ANI) custom tags directly to alignments table
        UPDATE alignments a
        SET 
            -- Add ZA tag with ANI value (ZA for "Z-ANI")
            tag_ZA = ROUND(t._ani, 4),
            
            -- Add ZS tag with shifted score (ZS for "Z-Score")
            tag_ZS = CASE 
                WHEN t.query_length > 0 THEN 
                    ROUND((t._alignment_score - (SELECT min_score FROM min_score) + 1.0) / CAST(t.query_length AS FLOAT), 4)
                ELSE NULL
            END
        FROM temp_alignment_calculations t
        WHERE a.qname = t.qname AND a.pos = t.pos AND a.rname = t.rname AND a.flag = t.flag
          AND t.query_length > 0;
        """

        con.execute(shifted_scores_query)

        # Drop the temporary calculation table now that we've used it
        con.execute("DROP TABLE IF EXISTS temp_alignment_calculations")

        # Log some statistics about the scores
        log.info("Collecting score statistics...")
        try:
            score_stats = con.execute(
                """
                SELECT 
                    MIN(tag_ZS) as min_shifted_score,
                    MAX(tag_ZS) as max_shifted_score,
                    AVG(tag_ZS) as avg_shifted_score,
                    COUNT(*) as total_alignments,
                    SUM(CASE WHEN tag_ZS IS NULL THEN 1 ELSE 0 END) as null_scores,
                    MIN(tag_ZA) as min_ani,
                    MAX(tag_ZA) as max_ani,
                    AVG(tag_ZA) as avg_ani
                FROM alignments
            """
            ).fetchone()

            log.info(
                f"Score statistics: min={score_stats[0]:.4f}, max={score_stats[1]:.4f}, avg={score_stats[2]:.4f}"
            )
            log.info(
                f"ANI statistics: min={score_stats[5]:.4f}, max={score_stats[6]:.4f}, avg={score_stats[7]:.4f}"
            )
            log.info(
                f"Total alignments: {score_stats[3]:,}, alignments with NULL scores: {score_stats[4]:,}"
            )
        except Exception as e:
            log.debug(f"Could not fetch score statistics: {e}")

        # --- BEGIN ADDED DEBUG LOGGING ---
        if log.isEnabledFor(logging.DEBUG):
            debug_alignments_query = """
            SELECT
                qname,
                rname,
                raw_tags,
                tag_NM,
                tag_AS,
                tag_MD
            FROM alignments
            WHERE raw_tags IS NOT NULL OR tag_NM IS NOT NULL OR tag_AS IS NOT NULL OR tag_MD IS NOT NULL
            LIMIT 5;
            """
            log.debug(
                "Sample of alignments table with parsed tags (where any tag might be present):"
            )
            try:
                sample_alignments = con.execute(debug_alignments_query).fetchall()
                if sample_alignments:
                    for row_idx, row in enumerate(sample_alignments):
                        log.debug(f"Alignments Sample Row {row_idx + 1}:")
                        log.debug(f"  QNAME: {row[0]}, RNAME: {row[1]}")
                        log.debug(f"  Raw Tags (in alignments table): {row[2]}")
                        log.debug(
                            f"  tag_NM: {row[3]}, tag_AS: {row[4]}, tag_MD: {row[5]}"
                        )
                else:
                    log.debug(
                        "No alignments found with non-NULL raw_tags or example parsed tags for sampling."
                    )
            except Exception as e_debug_aln:
                log.debug(
                    f"Could not fetch sample from alignments for debugging: {e_debug_aln}"
                )
        # --- END ADDED DEBUG LOGGING ---

    except Exception as e:
        log.error(f"Failed during data loading or initial processing: {e}")
        raise
    finally:
        # Reset progress bar setting to default
        try:
            con.execute("SET enable_progress_bar=false")
        except Exception:
            pass
        try:
            log.debug(f"Ensuring DuckDB threads are reset to {original_db_threads}.")
            con.execute(f"SET threads = {original_db_threads}")
        except Exception as thread_reset_e:
            log.warning(f"Could not reset DuckDB threads: {thread_reset_e}")
