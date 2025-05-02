import os
import logging
import subprocess
import shutil
import duckdb
import psutil
from pathlib import Path
from bam_filter.utils import check_tmp_dir_exists
from bam_filter.sam_utils_db import (
    load_sam_file,
    format_sam_output,
)
from bam_filter.db_manager import DatabaseManager
import tempfile
import json

log = logging.getLogger("my_logger")


def check_samtools():
    """Check if samtools is available in PATH"""
    if not shutil.which("samtools"):
        log.error(
            "samtools not found in PATH. Please install samtools and make sure it's in your PATH."
        )
        exit(1)


def convert_bam_to_parquet(input_file, output_path, args):
    """Convert BAM or SAM file to Parquet or DuckDB format via SAM"""
    tmp_dir = check_tmp_dir_exists(args.tmp_dir)
    compression = args.compression
    input_path = Path(input_file)
    input_ext = input_path.suffix.lower()
    output_format = args.output_format

    log.info(f"Converting {input_file} to {output_format.upper()} format")
    samtools_path = shutil.which("samtools")

    # Determine input file path early for better temp directory setting
    if input_ext == ".bam":
        check_samtools()
        log.info("Converting BAM to SAM for processing...")
        temp_sam_file = tempfile.NamedTemporaryFile(
            mode="w+t", suffix=".sam", delete=False
        )
        subprocess.run(
            [samtools_path, "view", "-h", input_file],
            stdout=temp_sam_file,
            check=True,
        )
        temp_sam_file.close()
        sam_file_path = temp_sam_file.name
        needs_cleanup = True
    elif input_ext == ".sam":
        # For SAM files, use the file directly
        sam_file_path = str(input_path.resolve())
        needs_cleanup = False
    else:
        log.error(f"Unexpected input file type for conversion: {input_ext}")
        exit(1)

    # Get input file directory for temp storage on same drive
    input_dir = os.path.dirname(os.path.abspath(sam_file_path))
    if not os.path.exists(input_dir) or not os.access(input_dir, os.W_OK):
        log.warning(
            f"Input directory {input_dir} not writable, using provided temp dir"
        )
        input_dir = tmp_dir.name

    try:
        # Create database manager with appropriate settings
        with DatabaseManager(
            temp_dir=input_dir,
            threads=args.threads,
            enable_progress=True,
            memory_limit=args.memory_limit,
            max_memory_pct=args.max_memory_pct,
            max_temp_size=args.max_temp_size,
        ) as db:
            # Load SAM data
            log.info("Loading SAM file into database...")
            load_sam_file(db.con, sam_file_path)

            # Count records for reporting
            header_count = db.execute("SELECT COUNT(*) FROM header").fetchone()[0]
            alignment_count = db.execute("SELECT COUNT(*) FROM alignments").fetchone()[
                0
            ]
            log.info(
                f"Loaded {header_count} header lines and {alignment_count} alignments"
            )

            # Write output based on format
            if output_format == "parquet":
                log.info("Writing Parquet output...")
                db.execute(
                    """
                    COPY alignments TO ? (FORMAT PARQUET, COMPRESSION ?);
                    COPY header TO ? (FORMAT PARQUET, COMPRESSION ?);
                    """,
                    [output_path, compression, f"{output_path}.header", compression],
                )
            elif output_format == "duckdb":
                log.info("Writing DuckDB output...")
                # Use the database manager's utility method
                db.create_persistent_database(
                    output_path, ["alignments", "header", "qname_index", "rname_index"]
                )

        # Clean up temporary file if created
        if needs_cleanup and os.path.exists(sam_file_path):
            os.unlink(sam_file_path)

    finally:
        tmp_dir.cleanup()

    log.info(f"Conversion complete. Output saved to {output_path}")


def convert_parquet_to_bam(input_parquet, output_bam, args):
    """Convert Parquet or DuckDB file to SAM format"""
    log.info(f"Converting {input_parquet} to SAM format")
    input_path = Path(input_parquet)
    input_format = "duckdb" if input_path.suffix.lower() == ".db" else "parquet"

    try:
        # Use database manager for in-memory database
        with DatabaseManager(
            enable_progress=True,
            memory_limit=args.memory_limit,
            max_memory_pct=args.max_memory_pct,
            max_temp_size=args.max_temp_size,
            threads=args.threads,
        ) as db:
            # Write SAM format directly
            output_sam_path = output_bam.replace(".bam", ".sam")
            log.info(f"Writing SAM header to {output_sam_path}")

            with open(output_sam_path, "w") as sam_out:
                # Read header based on input format
                if input_format == "parquet":
                    header_rows = db.execute(
                        f"SELECT header_line FROM read_parquet('{input_parquet}.header')"
                    ).fetchall()
                else:  # duckdb
                    source_db = duckdb.connect(database=str(input_path), read_only=True)
                    header_rows = source_db.execute(
                        "SELECT header_line FROM header"
                    ).fetchall()

                # Write header
                for row in header_rows:
                    sam_out.write(row[0] + "\n")

                # Get alignment data based on input format
                if input_format == "parquet":
                    alignment_query = format_sam_output(db.con)
                    result_relation = db.execute(
                        f"SELECT * FROM read_parquet('{input_parquet}')"
                    )
                else:  # duckdb
                    alignment_query = format_sam_output(source_db)
                    result_relation = source_db.execute("SELECT * FROM alignments")

                # Fetch and write in chunks for memory efficiency
                log.info("Streaming alignments to SAM file...")
                chunk_size = 100000
                processed_rows = 0
                while True:
                    chunk = result_relation.fetchmany(chunk_size)
                    if not chunk:
                        break
                    lines = ["\t".join(map(str, row)) + "\n" for row in chunk]
                    sam_out.writelines(lines)
                    processed_rows += len(chunk)

                log.info(f"Finished writing {processed_rows} alignment rows.")

                if input_format == "duckdb":
                    source_db.close()

    except Exception as e:
        log.error(f"Error converting file: {e}")
        raise

    log.info(f"Conversion complete. Output saved to {output_sam_path}")


def convert_files(args):
    """Convert between SAM/BAM and Parquet/DuckDB formats"""
    input_path = Path(args.input)
    input_ext = input_path.suffix.lower()

    # Supported input formats
    supported_input_formats = {".bam", ".sam", ".parquet", ".db"}
    if input_ext not in supported_input_formats:
        log.error(
            f"Unsupported input format: {input_ext}. Supported formats are: {', '.join(supported_input_formats)}"
        )
        exit(1)

    # Determine output path and extension
    if args.output:
        output_path = Path(args.output)
    else:
        # Auto-detect output format based on input and selected output format
        if input_ext in {".bam", ".sam"}:
            if args.output_format == "parquet":
                output_path = input_path.with_suffix(".parquet")
            else:  # duckdb (default)
                output_path = input_path.with_suffix(".db")
        elif input_ext in {".parquet", ".db"}:
            # Default output to SAM when converting from Parquet/DuckDB
            output_path = input_path.with_suffix(".sam")

    output_ext = output_path.suffix.lower()

    # Determine conversion direction
    if input_ext in {".bam", ".sam"}:
        if args.output_format == "parquet":
            log.info(
                f"Converting {input_ext.upper()[1:]} to Parquet: {input_path} -> {output_path}"
            )
        else:
            log.info(
                f"Converting {input_ext.upper()[1:]} to DuckDB: {input_path} -> {output_path}"
            )

        convert_bam_to_parquet(args.input, output_path, args)
    elif input_ext in {".parquet", ".db"}:
        if output_ext == ".bam":
            log.error(
                "Direct conversion to BAM is not supported. Convert to SAM first."
            )
            log.info(
                f"To convert to BAM, run: filterBAM convert -i {input_path} -o {output_path.with_suffix('.sam')}"
            )
            log.info(
                f"Then run: samtools view -b {output_path.with_suffix('.sam')} > {output_path}"
            )
            exit(1)
        elif output_ext == ".sam":
            source_format = "DuckDB" if input_ext == ".db" else "Parquet"
            log.info(
                f"Converting {source_format} to SAM: {input_path} -> {output_path}"
            )
            convert_parquet_to_bam(args.input, str(output_path), args)
        else:
            # Should not happen if logic above is correct, but as a safeguard
            log.error(f"Unsupported output format determined: {output_ext}")
            exit(1)
