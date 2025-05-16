import os
import logging
import duckdb
import psutil
from pathlib import Path
import tempfile

log = logging.getLogger("my_logger")


class DatabaseManager:
    """Manager for DuckDB connections with standardized configuration and utility methods."""

    def __init__(
        self,
        database=None,
        temp_dir=None,
        threads=None,
        memory_limit=None,
        max_memory_pct=80,
        enable_progress=True,
    ):
        """
        Initialize a new DatabaseManager.

        Args:
            database (str): Path to persistent database or ":memory:" for in-memory database
            temp_dir (str): Directory to store temporary files (None = use system temp)
            threads (int): Number of threads to use (None = auto-detect)
            memory_limit (str): Memory limit with units like "4GB" (None = auto-detect)
            max_memory_pct (int): Percentage of system memory to use when auto-detecting
            enable_progress (bool): Enable progress bar
        """
        # Handle temporary directory first
        if temp_dir is None:
            self.temp_dir = tempfile.gettempdir()
        else:
            self.temp_dir = temp_dir

        # Create the temp directory if it doesn't exist
        if not os.path.exists(self.temp_dir):
            try:
                os.makedirs(self.temp_dir, exist_ok=True)
                log.info(f"Created temporary directory: {self.temp_dir}")
            except Exception as e:
                log.warning(f"Could not create temp directory {self.temp_dir}: {e}")
                self.temp_dir = tempfile.gettempdir()
                log.info(f"Using system temp directory instead: {self.temp_dir}")

        # Check write permissions on temp directory
        if not os.access(self.temp_dir, os.W_OK):
            log.warning(f"No write permission on {self.temp_dir}")
            self.temp_dir = tempfile.gettempdir()
            log.info(f"Using system temp directory instead: {self.temp_dir}")

        # Always use a file-based database (never pure in‐memory)
        if database is None or database == ":memory:":
            # Create a temporary .db file in temp_dir
            # --- Add check and delete before creating ---
            try:
                tmp_db = tempfile.NamedTemporaryFile(
                    delete=False, suffix=".db", dir=self.temp_dir
                )
                tmp_db_path = tmp_db.name
                tmp_db.close()  # Close the file handle immediately

                # Now, ensure the file does not exist before DuckDB tries to use it
                if os.path.exists(tmp_db_path):
                    log.debug(
                        f"Deleting potentially pre-existing temp db file: {tmp_db_path}"
                    )
                    os.unlink(tmp_db_path)

                self.database = tmp_db_path  # Use the generated path
                self._temp_db_created = True

            except Exception as e:
                log.error(f"Failed to create or manage temporary database file: {e}")
                raise  # Re-raise the exception as this is critical

            # --- End check and delete ---

            if database == ":memory:":
                log.info(
                    f"Using temporary file-based database instead of in-memory: {self.database}"
                )
            else:
                log.info(
                    f"No database specified, created temporary file: {self.database}"
                )
        else:
            self.database = database
            self._temp_db_created = False
            log.info(f"Using specified database file: {self.database}")

        self.threads = threads if threads else max(4, os.cpu_count())

        # Auto-configure memory limit if not specified
        if memory_limit is None:
            available_mem = psutil.virtual_memory().available
            memory_bytes = int(available_mem * max_memory_pct / 100)
            memory_gb = max(1, memory_bytes // (1024 * 1024 * 1024))
            self.memory_limit = f"{memory_gb}GB"
            log.info(
                f"Auto-configured memory limit to {self.memory_limit} ({max_memory_pct}% of available memory)"
            )
        else:
            self.memory_limit = memory_limit
            log.info(f"Using specified memory limit: {self.memory_limit}")

        self.enable_progress = enable_progress
        self.con = None

    def __enter__(self):
        """Context manager entry point that connects to the database."""
        self.connect()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit point that closes the database connection."""
        if self.con:
            self.close()

        # Remove temporary DB file if we created it
        if getattr(self, "_temp_db_created", False):
            try:
                log.debug(f"Cleaning up temporary database file: {self.database}")
                os.unlink(self.database)
            except Exception as e:
                log.warning(
                    f"Could not remove temporary database file {self.database}: {e}"
                )

    def connect(self):
        """Establish and configure DuckDB connection."""
        if self.con:
            log.warning("Connection already established, reusing existing connection")
            return self.con

        log.debug(f"Connecting to DuckDB database: {self.database}")
        # Pass enable_external_access directly in the config
        # This is the primary way to enable it, must be done at connection time.
        config = {"enable_external_access": True}
        try:
            self.con = duckdb.connect(database=self.database, config=config)
            log.info("Successfully connected with external access enabled via config.")
        except Exception as connect_e:
            log.error(f"Failed to connect to database {self.database}: {connect_e}")
            # Try connecting without the config if it fails (older versions might not support it)
            try:
                log.warning(
                    "Retrying connection without external_access in config (external features may be limited)..."
                )
                self.con = duckdb.connect(database=self.database)
                log.warning(
                    "Connected without explicit external access. PROGRAM feature might not work."
                )
            except Exception as retry_e:
                log.error(f"Failed to connect to database on retry: {retry_e}")
                raise  # Re-raise the final connection error

        # Configure connection parameters
        self._configure_connection()
        return self.con

    def _configure_connection(self):
        """Apply configuration settings to the connection."""
        try:
            # Set critical settings first
            self.con.execute(f"SET threads TO {self.threads}")
            self.con.execute(f"SET memory_limit='{self.memory_limit}'")
            self.con.execute("PRAGMA enable_profiling")

            # Add optimizations for batch processing
            self.con.execute(
                "SET prefer_range_joins=true"
            )  # Prefer range joins when possible
            self.con.execute("SET explain_output='all'")  # Detailed execution plans
            self.con.execute(
                "SET force_compression='none'"
            )  # Disable compression during processing
            self.con.execute("SET preserve_insertion_order=false")
            self.con.execute(
                "SET immediate_transaction_mode=true"
            )  # Speed up transactions

            # Enable progress bar if requested
            if self.enable_progress:
                try:
                    self.con.execute(
                        "SET progress_bar_time=1000"
                    )  # Update every second
                    self.con.execute("SET enable_progress_bar=true")
                    log.info("Enabled progress bar for long-running queries")
                except Exception as e:
                    log.debug(f"Could not enable progress bar: {e}")

            # Set temp directory - critical for large dataset processing
            if self.temp_dir:
                log.info(f"Setting DuckDB temp directory to: {self.temp_dir}")
                safe_temp_dir = str(Path(self.temp_dir).resolve()).replace("'", "''")
                self.con.execute(f"SET temp_directory='{safe_temp_dir}'")

            # Apply performance optimizations that work across versions
            try:
                self.con.execute("SET checkpoint_threshold='1GB'")
                log.debug("Set checkpoint_threshold to 1GB")
            except Exception as e:
                log.debug(f"Could not set checkpoint_threshold: {e}")

            try:
                self.con.execute("SET allocator_flush_threshold='1GB'")
                log.debug("Set allocator_flush_threshold to 1GB")
            except Exception as e:
                log.debug(f"Could not set allocator_flush_threshold: {e}")

            try:
                self.con.execute("SET streaming_buffer_size='128MB'")
                log.debug("Set streaming_buffer_size to 128MB")
            except Exception as e:
                log.debug(f"Could not set streaming_buffer_size: {e}")

            # Try object cache which should work in most versions
            try:
                self.con.execute("SET enable_object_cache=true")
                log.debug("Enabled object cache")
            except Exception as e:
                log.debug(f"Could not enable object cache: {e}")

        except Exception as e:
            log.error(f"Error during basic DuckDB configuration: {e}")
            raise

    def close(self):
        """Close the database connection."""
        if self.con:
            self.con.close()
            self.con = None

    def execute(self, query, params=None):
        """Execute a query with optional parameters."""
        if not self.con:
            self.connect()

        try:
            if params:
                return self.con.execute(query, params)
            else:
                return self.con.execute(query)
        except Exception as e:
            log.error(f"Error executing query: {e}")
            log.debug(f"Query: {query}")
            if params:
                log.debug(f"Parameters: {params}")
            raise

    def inspect_table_storage(self, table_name, sample_rows=10):
        """
        Inspects the storage information for a given table using pragma_storage_info.

        Args:
            table_name (str): The name of the table to inspect.
            sample_rows (int): The number of rows to sample for the pragma.

        Returns:
            list: A list of tuples containing the storage information, or None if an error occurs.
        """
        if not self.con:
            self.connect()

        # Ensure table_name is properly quoted to prevent SQL injection if it comes from unsafe sources,
        # though typically table names are controlled. For pragma_storage_info, it needs to be a string literal.
        # DuckDB's Python API handles parameterization for values, but not for identifiers like table names in pragmas directly.
        # So, we construct the string carefully.
        if (
            not table_name.isalnum() and "_" not in table_name
        ):  # Basic check for valid table name characters
            log.error(f"Invalid table name for inspect_table_storage: {table_name}")
            return None

        query = f"""
        SELECT * EXCLUDE (column_path, segment_id, start, stats, persistent, block_id, block_offset, has_updates)
        FROM pragma_storage_info('{table_name}')
        USING SAMPLE {sample_rows} ROWS
        ORDER BY row_group_id;
        """
        log.debug(
            f"Executing storage inspection query for table '{table_name}':\n{query}"
        )
        try:
            result = self.con.execute(query).fetchall()
            return result
        except Exception as e:
            log.error(f"Error inspecting storage for table {table_name}: {e}")
            return None

    def inspect_all_managed_tables_storage(self, sample_rows=10):
        """
        Inspects and logs storage information for all standard managed tables.
        Standard tables are: 'alignments', 'reads', 'refs', 'header'.

        Args:
            sample_rows (int): The number of rows to sample for the pragma.
        """
        if not self.con:
            self.connect()

        managed_tables = ["alignments", "reads", "refs", "header"]
        log.info("--- Table Storage Inspection Start ---")
        for table_name in managed_tables:
            try:
                # Check if table exists before trying to inspect
                # No special quoting needed for 'refs' as it's not a keyword.
                query_table_name = table_name

                # Try to execute a query that would fail if the table doesn't exist
                # Using DESCRIBE is a common way to check for table existence and schema
                self.con.execute(f"DESCRIBE {query_table_name};").fetchall()
                # If DESCRIBE succeeds, the table exists.

                log.info(f"Inspecting storage for table: '{table_name}'")
                # The inspect_table_storage method itself passes table_name as a string literal to pragma_storage_info,
                # so no special quoting is needed for that call.
                storage_info = self.inspect_table_storage(table_name, sample_rows)
                if storage_info:
                    log.info(f"Storage information for '{table_name}':")
                    for row in storage_info:
                        log.info(f"  {row}")
                elif (
                    storage_info is None
                ):  # inspect_table_storage returned None due to an error within that method
                    log.warning(
                        f"Could not retrieve storage info for table '{table_name}' due to an error during its specific inspection."
                    )
                else:  # inspect_table_storage returned empty list
                    log.info(
                        f"No storage information returned by pragma for table '{table_name}'. It might be empty or not yet optimized for storage info."
                    )

            except Exception as e:
                # This catch is for errors like table not existing if DESCRIBE fails
                log.warning(f"Could not inspect storage for table '{table_name}': {e}")
        log.info("--- Table Storage Inspection End ---")

    def create_persistent_database(self, output_path, source_tables=None):
        """
        Create a persistent DuckDB database from the current in-memory database.

        Args:
            output_path (str): Path for the new persistent database
            source_tables (list): List of tables to copy (None = all tables)
        """
        if not self.con:
            raise ValueError("No active connection to copy from")

        # Remove the check for in-memory database since we always use file-based now
        log.info(f"Creating persistent database at {output_path}")
        # Ensure output path is absolute for ATTACH
        abs_output_path = str(Path(output_path).resolve())
        # Ensure source path is absolute for ATTACH
        abs_source_path = str(Path(self.database).resolve())

        # --- Delete existing target database file before creating ---
        if os.path.exists(abs_output_path):
            try:
                os.unlink(abs_output_path)
                log.debug(f"Deleted existing target database file: {abs_output_path}")
            except OSError as e:
                log.error(
                    f"Could not delete existing database file {abs_output_path}: {e}"
                )
                raise  # Re-raise error if deletion fails
        # --- End deletion ---

        # Connect to the target database (or create it)
        out_db = duckdb.connect(database=abs_output_path)

        # Configure the output database for optimal writing performance - using try/except for each setting
        try:
            out_db.execute(f"SET threads = {self.threads}")
        except Exception as e:
            log.debug(f"Could not set threads for output db: {e}")

        try:
            out_db.execute("SET preserve_insertion_order=false")
        except Exception as e:
            log.debug(f"Could not set preserve_insertion_order for output db: {e}")

        try:
            out_db.execute("SET enable_object_cache=true")
        except Exception as e:
            log.debug(f"Could not set enable_object_cache for output db: {e}")

        try:
            out_db.execute("SET checkpoint_on_shutdown=true")
        except Exception as e:
            log.debug(f"Could not set checkpoint_on_shutdown for output db: {e}")

        source_db_alias = "source_db"
        try:
            # Attach the source database (current connection's file) to the target connection
            log.debug(
                f"Attaching source database '{abs_source_path}' as '{source_db_alias}'"
            )
            out_db.execute(
                f"ATTACH '{abs_source_path}' AS {source_db_alias} (READ_ONLY)"
            )

            # If no specific tables are provided, get all tables from the source database
            if source_tables is None:
                tables_result = out_db.execute(
                    f"SELECT name FROM {source_db_alias}.sqlite_master WHERE type='table'"
                ).fetchall()
                # Filter out internal/sqlite tables if necessary, though DuckDB usually handles this
                source_tables = [
                    row[0] for row in tables_result if not row[0].startswith("sqlite_")
                ]
                log.debug(f"Auto-detected tables to copy: {source_tables}")

            # Copy each table and create indexes
            for table in source_tables:
                log.info(f"Copying table {table} to persistent database")
                quoted_target_table = f'"{table}"'

                # Use faster COPY approach for large tables
                if table in ["alignments", "alignment_references"]:
                    # For large tables, use a more optimized approach
                    out_db.execute(
                        f"CREATE TABLE {quoted_target_table} AS SELECT * FROM {source_db_alias}.{table} LIMIT 0"
                    )
                    out_db.execute(
                        f"INSERT INTO {quoted_target_table} SELECT * FROM {source_db_alias}.{table}"
                    )
                else:
                    # For smaller tables, use the simpler CREATE TABLE AS approach
                    out_db.execute(
                        f"CREATE TABLE {quoted_target_table} AS SELECT * FROM {source_db_alias}.{table}"
                    )

                # Create typical indexes based on table name patterns
                self._create_indexes_for_table(out_db, table)

                # Force intermediate checkpoint for large tables
                if table in ["alignments", "alignment_references"]:
                    out_db.execute("CHECKPOINT")

        except Exception as e:
            log.error(f"Error during persistent database creation: {e}")
            raise
        finally:
            # Detach the source database and close the target connection
            try:
                log.debug(f"Detaching source database '{source_db_alias}'")
                out_db.execute(f"DETACH {source_db_alias}")
            except Exception as detach_e:
                log.warning(
                    f"Could not detach database '{source_db_alias}': {detach_e}"
                )
            finally:
                out_db.close()

        log.info(f"Persistent database created at {abs_output_path}")
        return abs_output_path

    def _create_indexes_for_table(self, con_target, table_name):
        """Internal helper to create standard indexes for a given table."""
        # Most indexes are now unnecessary as DuckDB can auto-optimize
        # Only create special-case indexes for specific query patterns
        log.debug(f"Analyzing table {table_name} for query optimization")
        try:
            quoted_table = f'"{table_name}"'  # Ensure table name is quoted for safety, though 'refs' doesn't strictly need it.

            # Run ANALYZE to collect statistics for better query planning
            con_target.execute(f"ANALYZE {quoted_table}")

            # Only create special indexes for very specific access patterns
            # that DuckDB might not optimize automatically
            if table_name == "refs":
                # The ref_name column is frequently used for joining/filtering
                con_target.execute(
                    f"CREATE INDEX IF NOT EXISTS idx_ref_name ON {quoted_table}(ref_name)"
                )
        except Exception as index_e:
            log.warning(f"Could not analyze table {table_name}: {index_e}")

    def add_temp_management_methods(self):
        """Add a method to explain temp file usage"""
        temp_usage_info = """
DuckDB Temporary File Usage Information:

1. Spilling Data: When processing large datasets exceeds available memory, DuckDB "spills" 
   intermediate results to disk to continue processing.
   
2. Join Operations: Complex joins on large tables often require materialization of temporary data.

3. Group By & Aggregates: When grouping large datasets, temporary space is used for hash tables.

4. Window Functions: Can generate large sorted temporary results on disk.

5. Order By: Sorting large datasets requires temporary storage.

Optimization Tips:
- Increase memory_limit if RAM is available
- Process data in smaller chunks
- Use simpler queries where possible
- Apply filtering early in the query to reduce intermediate data size
- Use temp_directory_compression to compress temporary files
- Set max_temp_directory_size to limit disk usage (may cause queries to fail)
"""
        return temp_usage_info

    def set_max_memory(self, memory_gb):
        """Update memory limit during runtime"""
        if self.con:
            try:
                self.con.execute(f"SET memory_limit='{memory_gb}GB'")
                self.memory_limit = f"{memory_gb}GB"
                log.info(f"Updated memory limit to {memory_gb}GB")
                return True
            except Exception as e:
                log.error(f"Failed to update memory limit: {e}")
                return False
        return False

    def export_temp_filtered_header(self, output_path):
        """Export the temp_filtered_header table to a Parquet file."""
        try:
            self.con.execute(
                f"COPY temp_filtered_header TO '{output_path}/header.parquet' "
                f"(FORMAT PARQUET, COMPRESSION zstd, ROW_GROUP_SIZE 100000)"
            )
            log.info(f"Exported temp_filtered_header to {output_path}/header.parquet")
        except Exception as e:
            log.error(f"Failed to export temp_filtered_header: {e}")
            raise

    def export_optimized_table(self, table_name, output_path, compression_level=11):
        """
        Export a table with optimized settings based on its type.

        Args:
            table_name (str): Table to export ('alignments', 'reads', 'refs', 'header')
            output_path (str): Directory to save the Parquet file
            compression_level (int): zstd compression level (1-22, higher = smaller but slower)

        Returns:
            bool: True if export successful, False otherwise
        """
        from bam_filter.sam_utils_db import export_optimized_table

        return export_optimized_table(
            self.con, table_name, output_path, compression_level
        )

    def get_top_alignments(
        self, min_ani=0.9, min_read_length=50, max_read_length=None, limit=1000
    ):
        """
        Get top alignments based on shifted score.

        Args:
            min_ani (float): Minimum ANI (Average Nucleotide Identity) threshold (0-1)
            min_read_length (int): Minimum read length to consider
            max_read_length (int, optional): Maximum read length to consider
            limit (int): Maximum number of results to return

        Returns:
            List[Dict]: List of alignment dictionaries with top scores
        """
        if not self.con:
            self.connect()

        query_conditions = [f"ani >= {min_ani}", f"query_length >= {min_read_length}"]

        if max_read_length:
            query_conditions.append(f"query_length <= {max_read_length}")

        where_clause = " AND ".join(query_conditions)

        query = f"""
        WITH top_alignments AS (
            SELECT
                qname,
                rname,
                query_length,
                alignment_score,
                shifted_score,
                ani,
                ROW_NUMBER() OVER (PARTITION BY qname ORDER BY shifted_score DESC) as rank
            FROM
                alignments
            WHERE
                {where_clause}
        )
        SELECT
            qname,
            rname,
            query_length,
            alignment_score,
            shifted_score,
            ani
        FROM
            top_alignments
        WHERE
            rank = 1
        ORDER BY
            shifted_score DESC
        LIMIT {limit};
        """

        try:
            results = self.con.execute(query).fetchall()
            return [
                {
                    "query_id": row[0],
                    "subject_id": row[1],
                    "query_length": row[2],
                    "alignment_score": row[3],
                    "shifted_score": row[4],
                    "ani": row[5],
                }
                for row in results
            ]
        except Exception as e:
            log.error(f"Error retrieving top alignments: {e}")
            return []

    def cleanup_temp_files(self):
        """Clean up temporary files that might have been created during database operations."""
        log.debug("Cleaning up any temporary files from database operations...")
        try:
            # Remove checkpoint call, just check connection
            if self.con:
                log.debug("Connected database ready for cleanup")

            # Set temp_directory may have created temp files that need special cleanup
            if hasattr(self, "temp_dir") and self.temp_dir:
                log.debug(f"Using temp dir for cleanup check: {self.temp_dir}")
                # We can't directly delete temp files as they may still be in use
                # Just log that we're ensuring proper release of resources

            return True
        except Exception as e:
            log.warning(f"Error during temp file cleanup: {e}")
            return False

    def export_optimized_alignments(self, output_dir_str, compression_level=11):
        """Export alignments table with optimized settings specifically for it."""
        from .sam_utils_db import export_optimized_table

        return export_optimized_table(
            self.con, "alignments", output_dir_str, compression_level
        )

    def table_exists(self, table_name):
        """Check if a table exists in the database.

        Args:
            table_name (str): Name of the table to check

        Returns:
            bool: True if the table exists, False otherwise
        """
        try:
            result = self.con.execute(
                f"SELECT COUNT(1) FROM information_schema.tables WHERE table_name='{table_name}'"
            ).fetchone()
            return result[0] > 0
        except Exception as e:
            self.log.debug(f"Error checking if table '{table_name}' exists: {e}")
            # Alternative approach if the above fails
            try:
                self.con.execute(f"SELECT * FROM {table_name} LIMIT 0")
                return True
            except Exception:
                return False

    def get_table_columns(self, table_name):
        """Get the column names for a table.

        Args:
            table_name (str): Name of the table

        Returns:
            list: List of column names
        """
        try:
            columns = self.con.execute(f"PRAGMA table_info('{table_name}')").fetchall()
            return [col[1] for col in columns]
        except Exception as e:
            self.log.debug(f"Error getting columns for table '{table_name}': {e}")
            return []
