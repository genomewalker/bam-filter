import logging
import numpy as np
from typing import Dict, NamedTuple, Optional
import psutil

log = logging.getLogger("my_logger")


class ColumnInfo(NamedTuple):
    """Information about a column's data type and properties"""

    numpy_type: str
    nullable: bool = True


# Define mapping from DuckDB types to NumPy types
DUCKDB_TO_NUMPY_TYPES = {
    "SMALLINT": "int16",
    "INTEGER": "int32",
    "BIGINT": "int64",
    "REAL": "float32",
    "DOUBLE": "float64",
    "VARCHAR": "object",  # For strings
    "BOOLEAN": "bool",
}


def get_column_info(duckdb_type: str) -> ColumnInfo:
    """Convert DuckDB type to ColumnInfo with appropriate NumPy type"""
    # Remove any length specifiers and normalize
    base_type = duckdb_type.split("(")[0].upper()
    numpy_type = DUCKDB_TO_NUMPY_TYPES.get(base_type, "object")
    return ColumnInfo(numpy_type=numpy_type)


def get_available_memory() -> int:
    """Get available system memory in bytes"""
    try:
        return psutil.virtual_memory().available
    except:
        # Fallback to a conservative estimate if psutil fails
        return 1024 * 1024 * 1024  # 1GB


def calculate_optimal_row_group_size(
    con, table_name, columns_info_list, available_memory: Optional[int] = None
) -> int:
    """
    Calculate optimal row group size based on Parquet recommendations (512MB-1GB per group)

    Args:
        con: DuckDB connection object
        table_name: Name of the table
        columns_info_list: List of ColumnInfo objects
        available_memory: Optional memory limit in bytes

    Returns:
        int: Optimal row group size
    """
    if available_memory is None:
        available_memory = get_available_memory()

    # Get the total rows in the table
    total_rows = con.execute(f'SELECT COUNT(*) FROM "{table_name}"').fetchone()[0]

    # Calculate average row size based on column types
    row_size = 0

    # Handle columns_info_list whether it's a list of ColumnInfo objects or a dict
    if isinstance(columns_info_list, list):
        # Process as a list of ColumnInfo objects
        for col_info in columns_info_list:
            # Get size in bytes for each type
            if hasattr(col_info, "numpy_type"):
                dtype = np.dtype(col_info.numpy_type)
                row_size += dtype.itemsize
            elif isinstance(col_info, str):
                # Handle if we somehow got a list of strings
                log.warning(
                    f"Received string '{col_info}' in columns_info_list instead of ColumnInfo object. Using default size."
                )
                row_size += 8  # Use a reasonable default (8 bytes)
    elif isinstance(columns_info_list, dict):
        # Process as a dict of ColumnInfo objects (original behavior)
        for col_info in columns_info_list.values():
            dtype = np.dtype(col_info.numpy_type)
            row_size += dtype.itemsize
    else:
        # Fallback if we got something unexpected
        log.warning(
            f"Unexpected columns_info_list type: {type(columns_info_list)}. Using default row size."
        )
        row_size = 16  # Use a reasonable default for a single column

    # Ensure we have a non-zero row size
    if row_size <= 0:
        log.warning(f"Calculated row_size was {row_size}, using default of 16 bytes.")
        row_size = 16  # Use a reasonable default size

    # Target row group size (aim for 768MB = middle of recommended range)
    TARGET_GROUP_SIZE = 768 * 1024 * 1024  # 768MB in bytes

    # Calculate rows needed to reach target size
    rows_for_target = TARGET_GROUP_SIZE // row_size

    # Set bounds based on Parquet recommendations
    MIN_GROUP_SIZE_BYTES = 512 * 1024 * 1024  # 512MB
    MAX_GROUP_SIZE_BYTES = 1024 * 1024 * 1024  # 1GB

    min_rows = MIN_GROUP_SIZE_BYTES // row_size
    max_rows = MAX_GROUP_SIZE_BYTES // row_size

    # Ensure we don't exceed total rows
    max_rows = min(max_rows, total_rows) if total_rows > 0 else max_rows

    # Adjust based on available memory (ensure group can fit in memory)
    mem_limited_rows = (
        available_memory // 2
    ) // row_size  # Use at most 50% of available memory

    # Select final row count
    row_group_size = min(rows_for_target, mem_limited_rows, max_rows)
    row_group_size = max(row_group_size, min_rows)

    # Log the decision
    group_size_mb = (row_group_size * row_size) / (1024 * 1024)
    log.debug(
        f"""
        Row group size calculation for table '{table_name}':
        - Row size: {row_size} bytes
        - Total rows: {total_rows:,}
        - Target group size: 768 MB
        - Calculated rows per group: {row_group_size:,}
        - Actual group size: {group_size_mb:.2f} MB
        - Total groups: {total_rows / max(1, row_group_size):.1f}
    """
    )

    return int(row_group_size)
