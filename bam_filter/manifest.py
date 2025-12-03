#!/usr/bin/env python3
"""
Provenance tracking for Parquet conversions and operations.

This module provides functionality to track all operations performed on
Parquet datasets, creating an append-only audit trail of transformations.
"""

import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Any, Optional, List


MANIFEST_VERSION = "1.0"
MANIFEST_FILENAME = "manifest.json"


def get_filterbam_version() -> str:
    """Get the current filterBAM version."""
    try:
        # Try to get version from package
        from bam_filter import __version__
        return __version__
    except (ImportError, AttributeError):
        # Fallback to reading from git or setup.py
        try:
            import subprocess
            result = subprocess.run(
                ['git', 'describe', '--tags', '--always'],
                capture_output=True,
                text=True,
                check=True,
                cwd=Path(__file__).parent.parent
            )
            return result.stdout.strip()
        except:
            return "unknown"


def get_current_timestamp() -> str:
    """Get current UTC timestamp in ISO format."""
    return datetime.now(timezone.utc).isoformat()


def get_file_info(file_path: str) -> Dict[str, Any]:
    """Get information about a file."""
    path = Path(file_path)
    if not path.exists():
        return {
            "path": str(file_path),
            "exists": False
        }

    stat = path.stat()
    return {
        "path": str(file_path),
        "exists": True,
        "size_bytes": stat.st_size,
        "modified": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat()
    }


def create_manifest(
    output_dir: str,
    operation: str,
    command: str,
    parameters: Dict[str, Any],
    input_info: Dict[str, Any],
    output_info: Dict[str, Any],
    original_input: Optional[Dict[str, Any]] = None
) -> None:
    """
    Create a new manifest file with the initial operation.

    Args:
        output_dir: Directory where manifest will be saved
        operation: Operation name (e.g., "bam_to_parquet")
        command: Full command that was executed
        parameters: Operation parameters
        input_info: Information about input data
        output_info: Information about output/results
        original_input: Information about the original source file (optional)
    """
    manifest_path = Path(output_dir) / MANIFEST_FILENAME

    # Use original_input if provided, otherwise use input_info
    if original_input is None:
        original_input = input_info.copy()

    manifest = {
        "schema_version": MANIFEST_VERSION,
        "created": get_current_timestamp(),
        "original_input": original_input,
        "operations": [
            {
                "operation_id": 1,
                "timestamp": get_current_timestamp(),
                "operation": operation,
                "command": command,
                "filterbam_version": get_filterbam_version(),
                "parameters": parameters,
                "input": input_info,
                "output": output_info
            }
        ],
        "current_state": {
            "total_records": output_info.get("records_processed", 0),
            "last_modified": get_current_timestamp(),
            "operations_count": 1
        }
    }

    # Write manifest
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    print(f"Manifest created: {manifest_path}")


def append_operation(
    output_dir: str,
    operation: str,
    command: str,
    parameters: Dict[str, Any],
    input_info: Dict[str, Any],
    output_info: Dict[str, Any]
) -> None:
    """
    Append a new operation to an existing manifest.

    Args:
        output_dir: Directory containing the manifest
        operation: Operation name (e.g., "filter", "reassign")
        command: Full command that was executed
        parameters: Operation parameters
        input_info: Information about input data
        output_info: Information about output/results
    """
    manifest_path = Path(output_dir) / MANIFEST_FILENAME

    if not manifest_path.exists():
        raise FileNotFoundError(
            f"No manifest found at {manifest_path}. "
            "Cannot append operation to non-existent manifest."
        )

    # Read existing manifest
    with open(manifest_path, 'r') as f:
        manifest = json.load(f)

    # Create new operation entry
    operation_id = len(manifest["operations"]) + 1
    new_operation = {
        "operation_id": operation_id,
        "timestamp": get_current_timestamp(),
        "operation": operation,
        "command": command,
        "filterbam_version": get_filterbam_version(),
        "parameters": parameters,
        "input": input_info,
        "output": output_info
    }

    # Append operation
    manifest["operations"].append(new_operation)

    # Update current state
    manifest["current_state"]["total_records"] = output_info.get(
        "records_processed",
        manifest["current_state"]["total_records"]
    )
    manifest["current_state"]["last_modified"] = get_current_timestamp()
    manifest["current_state"]["operations_count"] = operation_id

    # Write updated manifest
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    print(f"Manifest updated: {manifest_path} (operation {operation_id}: {operation})")


def read_manifest(output_dir: str) -> Dict[str, Any]:
    """
    Read and return the manifest.

    Args:
        output_dir: Directory containing the manifest

    Returns:
        Manifest dictionary

    Raises:
        FileNotFoundError: If manifest doesn't exist
    """
    manifest_path = Path(output_dir) / MANIFEST_FILENAME

    if not manifest_path.exists():
        raise FileNotFoundError(f"No manifest found at {manifest_path}")

    with open(manifest_path, 'r') as f:
        return json.load(f)


def validate_manifest(output_dir: str) -> bool:
    """
    Validate manifest structure and consistency.

    Args:
        output_dir: Directory containing the manifest

    Returns:
        True if valid, False otherwise
    """
    try:
        manifest = read_manifest(output_dir)

        # Check required fields
        required_fields = ["schema_version", "created", "original_input", "operations", "current_state"]
        for field in required_fields:
            if field not in manifest:
                print(f"Manifest validation failed: missing field '{field}'")
                return False

        # Check operations are sequential
        operations = manifest["operations"]
        for i, op in enumerate(operations, 1):
            if op["operation_id"] != i:
                print(f"Manifest validation failed: operation IDs not sequential")
                return False

        # Check current state matches operations
        if manifest["current_state"]["operations_count"] != len(operations):
            print(f"Manifest validation failed: operations count mismatch")
            return False

        return True

    except Exception as e:
        print(f"Manifest validation failed: {e}")
        return False


def print_manifest_summary(output_dir: str) -> None:
    """
    Print a human-readable summary of the manifest.

    Args:
        output_dir: Directory containing the manifest
    """
    try:
        manifest = read_manifest(output_dir)

        print("="*80)
        print("MANIFEST SUMMARY")
        print("="*80)
        print(f"Created: {manifest['created']}")
        print(f"Original input: {manifest['original_input'].get('path', 'N/A')}")
        print(f"Operations: {manifest['current_state']['operations_count']}")
        print(f"Current records: {manifest['current_state']['total_records']:,}")
        print(f"Last modified: {manifest['current_state']['last_modified']}")
        print()

        print("OPERATION HISTORY:")
        print("-"*80)
        for op in manifest["operations"]:
            print(f"{op['operation_id']}. {op['operation']} ({op['timestamp']})")
            print(f"   Command: {op['command']}")
            if 'records_processed' in op['output']:
                print(f"   Records: {op['output']['records_processed']:,}")
            print()

        print("="*80)

    except Exception as e:
        print(f"Failed to read manifest: {e}")


def get_command_line() -> str:
    """
    Reconstruct the command line that was executed.

    Returns:
        Command line string
    """
    # Get the command line arguments
    return ' '.join(sys.argv)


if __name__ == '__main__':
    # Test/demo code
    import argparse

    parser = argparse.ArgumentParser(description="Manifest utilities")
    parser.add_argument('action', choices=['validate', 'summary'], help="Action to perform")
    parser.add_argument('output_dir', help="Output directory containing manifest")

    args = parser.parse_args()

    if args.action == 'validate':
        if validate_manifest(args.output_dir):
            print("✓ Manifest is valid")
            sys.exit(0)
        else:
            print("✗ Manifest is invalid")
            sys.exit(1)

    elif args.action == 'summary':
        print_manifest_summary(args.output_dir)
