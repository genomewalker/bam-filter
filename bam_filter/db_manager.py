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
        database=":memory:",
        temp_dir=None,
        threads=None,
        memory_limit=None,
        max_memory_pct=80,
        enable_progress=True,
        max_temp_size="500GB",
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
            max_temp_size (str): Maximum temporary directory size (with units like "500GB")
        """
        self.database = database

        # Handle temporary directory
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
        self.max_temp_size = max_temp_size
        self.con = None

    def __enter__(self):
        """Context manager entry point that connects to the database."""
        self.connect()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit point that closes the database connection."""
        if self.con:
            self.close()

    def connect(self):
        """Establish and configure DuckDB connection."""
        if self.con:
            log.warning("Connection already established, reusing existing connection")
            return self.con

        log.debug(f"Connecting to DuckDB database: {self.database}")
        self.con = duckdb.connect(database=self.database)

        # Configure connection parameters
        self._configure_connection()
        return self.con

    def _configure_connection(self):
        """Apply configuration settings to the connection."""
        try:
            if self.enable_progress:
                self.con.execute("PRAGMA enable_progress_bar=true")

            # Set temp directory - critical for large dataset processing
            if self.temp_dir:
                log.info(f"Setting DuckDB temp directory to: {self.temp_dir}")
                self.con.execute(f"SET temp_directory='{self.temp_dir}'")

            # Set thread count
            log.info(f"Setting DuckDB thread count to: {self.threads}")
            self.con.execute(f"SET threads TO {self.threads}")

            # Set memory limit
            log.info(f"Setting DuckDB memory limit to: {self.memory_limit}")
            self.con.execute(f"SET memory_limit='{self.memory_limit}'")

            # Set maximum temp directory size
            log.info(f"Setting DuckDB max temp directory size to: {self.max_temp_size}")
            self.con.execute(f"PRAGMA max_temp_directory_size='{self.max_temp_size}'")

            # Add temp file specific optimizations
            try:
                # Enable compression for temporary data when spilling to disk
                self.con.execute("PRAGMA temp_directory_compression='zstd'")

                # Try to minimize disk usage by doing more work in memory
                self.con.execute("PRAGMA memory_limit_affinity='disk_to_memory'")

                # Clean up temp files more aggressively
                self.con.execute("PRAGMA cleanup_on_close=true")

                # Reduce materialization of intermediate results where possible
                self.con.execute("PRAGMA enable_intermediate_materialization=false")

                # Use compressed execution where possible
                self.con.execute("PRAGMA optimize_compressed_materialization=true")

                log.info(
                    "Configured temp storage with optimizations for reduced disk usage"
                )
            except Exception as e:
                log.debug(f"Some temp storage optimizations not supported: {e}")

            # Optimize for performance
            try:
                self.con.execute("PRAGMA preserve_insertion_order=false")
                self.con.execute("PRAGMA enable_object_cache")
                self.con.execute("PRAGMA enable_profiling")
                self.con.execute("PRAGMA memory_limit_affinity='memory_to_disk'")
            except Exception as e:
                log.debug(
                    f"Some optimizations not supported in this DuckDB version: {e}"
                )

        except Exception as e:
            log.error(f"Error configuring DuckDB connection: {e}")
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

    def create_persistent_database(self, output_path, source_tables=None):
        """
        Create a persistent DuckDB database from the current in-memory database.

        Args:
            output_path (str): Path for the new persistent database
            source_tables (list): List of tables to copy (None = all tables)
        """
        if not self.con:
            raise ValueError("No active connection to copy from")

        log.info(f"Creating persistent database at {output_path}")
        out_db = duckdb.connect(database=str(output_path))

        # If no specific tables are provided, get all tables from the current connection
        if source_tables is None:
            tables_result = self.con.execute("SHOW TABLES").fetchall()
            source_tables = [row[0] for row in tables_result]

        # Copy each table and create indexes
        for table in source_tables:
            log.info(f"Copying table {table} to persistent database")
            out_db.execute(f"CREATE TABLE {table} AS SELECT * FROM con.{table}")

            # Create typical indexes based on table name patterns
            if table == "alignments":
                out_db.execute(f"CREATE INDEX idx_{table}_qname ON {table}(qname)")
                out_db.execute(f"CREATE INDEX idx_{table}_rname ON {table}(rname)")
            elif table == "qname_index":
                out_db.execute(f"CREATE INDEX idx_qname ON {table}(qname)")
            elif table == "rname_index":
                out_db.execute(f"CREATE INDEX idx_rname ON {table}(rname)")

        out_db.close()
        log.info(f"Persistent database created at {output_path}")

        return output_path

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

    def cleanup_temp_files(self):
        """Force cleanup of temporary files"""
        if self.con:
            try:
                self.con.execute("PRAGMA force_cleanup()")
                log.info("Forced cleanup of temporary DuckDB files")
            except:
                log.warning(
                    "Could not force cleanup - may not be supported in this version"
                )


# Example usage patterns:

# Single query execution:
# db = DatabaseManager(temp_dir="/tmp")
# result = db.execute("SELECT 1")
# db.close()

# Multiple operations with context manager:
# with DatabaseManager(threads=4) as db:
#     db.execute("CREATE TABLE test (id INTEGER, name VARCHAR)")
#     db.execute("INSERT INTO test VALUES (1, 'test')")
#     result = db.execute("SELECT * FROM test").fetchall()

# Creating a persistent database:
# with DatabaseManager() as db:
#     db.execute("CREATE TABLE data (id INTEGER)")
#     db.execute("INSERT INTO data SELECT * FROM range(1000)")
#     db.create_persistent_database("output.db", ["data"])
