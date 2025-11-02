"""
High-level Python interface for ultra-fast taxonomy database operations.

This module provides a simple, user-friendly interface to the high-performance
Cython taxonomy database implementation.

Example usage:
    # Build database from NCBI dump files
    >>> from bam_filter.taxonomy import TaxonomyDB
    >>> tax = TaxonomyDB.from_ncbi('nodes.dmp', 'names.dmp', 'acc2taxid.txt')
    >>> tax.save('taxonomy_db')

    # Load pre-built database (100x faster)
    >>> tax = TaxonomyDB.load('taxonomy_db')

    # Query operations
    >>> tax.get_name(9606)
    'Homo sapiens'

    >>> tax.get_lineage(9606)
    [1, 131567, 2759, 33154, ..., 9605, 9606]

    >>> tax.compute_lca([9606, 9598])  # Human and Chimp
    9604  # Hominidae

    # Get taxid from accession
    >>> tax.acc2taxid('NC_000001.11')
    9606

    # Batch operations (much faster)
    >>> references = ['NC_000001.11', 'NC_000002.12', ...]
    >>> taxids = tax.acc2taxid_batch(references)
    >>> lca = tax.compute_lca(taxids)
"""

import os
from pathlib import Path
from typing import List, Optional, Union
import logging

from bam_filter.taxonomy_db import (
    TaxonomyDatabase,
    AccessionMapping,
    load_taxonomy_from_ncbi,
    load_accession_map_from_file,
)

log = logging.getLogger("my_logger")


class TaxonomyDB:
    """
    High-level interface for taxonomy database operations.

    This class combines taxonomy tree structure and accession mappings
    into a single convenient interface.
    """

    def __init__(self, tax_db: TaxonomyDatabase, acc_map: Optional[AccessionMapping] = None):
        """
        Initialize TaxonomyDB.

        Parameters
        ----------
        tax_db : TaxonomyDatabase
            Taxonomy database object
        acc_map : AccessionMapping, optional
            Accession to taxid mapping
        """
        self.tax_db = tax_db
        self.acc_map = acc_map

    @classmethod
    def from_ncbi(
        cls,
        nodes_file: str,
        names_file: str,
        acc2taxid_file: Optional[str] = None,
        is_custom: bool = False,
        num_threads: int = 1,
        build_lca_cache: bool = False,
        cache_taxids: Optional[List[int]] = None,
    ) -> 'TaxonomyDB':
        """
        Build taxonomy database from NCBI dump files.

        Parameters
        ----------
        nodes_file : str
            Path to nodes.dmp file
        names_file : str
            Path to names.dmp file
        acc2taxid_file : str, optional
            Path to acc2taxid file
        is_custom : bool, optional
            Whether using custom taxonomy format (default: False)
        num_threads : int, optional
            Number of threads for construction (default: 1)
        build_lca_cache : bool, optional
            Whether to build LCA cache for faster queries (default: False)
        cache_taxids : list of int, optional
            Specific taxids to cache (if None and build_lca_cache=True, caches all)

        Returns
        -------
        TaxonomyDB
            Initialized taxonomy database
        """
        log.info("Building taxonomy database from NCBI files...")

        # Load taxonomy tree
        tax_db = load_taxonomy_from_ncbi(nodes_file, names_file, num_threads=num_threads)

        # Load accession mapping if provided
        acc_map = None
        if acc2taxid_file is not None:
            acc_map = load_accession_map_from_file(
                acc2taxid_file,
                is_custom=is_custom,
                num_threads=num_threads
            )

        # Build LCA cache if requested
        if build_lca_cache:
            if cache_taxids is None:
                log.warning("LCA cache requested but no taxids specified. Cache will not be built.")
            else:
                log.info(f"Building LCA cache for {len(cache_taxids):,} taxids...")
                tax_db.build_lca_cache(cache_taxids)

        return cls(tax_db, acc_map)

    def load_additional_accessions(self, acc2taxid_file: str, num_threads: int = 1, is_custom: bool = False):
        """
        Load additional accession mappings from another acc2taxid file.

        This allows merging multiple acc2taxid files (e.g., nucl_gb + nucl_wgs + prot)
        into a single taxonomy database.

        Parameters
        ----------
        acc2taxid_file : str
            Path to additional acc2taxid file
        num_threads : int, optional
            Number of threads for parsing (default: 1)
        is_custom : bool, optional
            Whether using custom format (default: False)
        """
        import pandas as pd
        import tempfile
        from bam_filter.taxonomy_db import _build_accession_map_from_dataframe

        if self.acc_map is None:
            # First accession file - create new mapping
            log.info(f"Loading accessions from {acc2taxid_file}...")
            self.acc_map = load_accession_map_from_file(
                acc2taxid_file,
                is_custom=is_custom,
                num_threads=num_threads
            )
        else:
            # Merge with existing mappings using temporary Parquet files
            log.info(f"Merging accessions from {acc2taxid_file}...")

            # Load additional accessions
            additional_map = load_accession_map_from_file(
                acc2taxid_file,
                is_custom=is_custom,
                num_threads=num_threads
            )

            # Use temporary files to extract data
            with tempfile.TemporaryDirectory() as tmpdir:
                existing_parquet = os.path.join(tmpdir, 'existing.parquet')
                additional_parquet = os.path.join(tmpdir, 'additional.parquet')

                # Save both to Parquet
                self.acc_map.to_parquet(existing_parquet)
                additional_map.to_parquet(additional_parquet)

                # Load as DataFrames
                existing_df = pd.read_parquet(existing_parquet)
                additional_df = pd.read_parquet(additional_parquet)

                # Combine and deduplicate
                combined_df = pd.concat([existing_df, additional_df], ignore_index=True)
                combined_df = combined_df.drop_duplicates(subset='accession', keep='first')

                log.info(f"Merged {len(additional_df):,} accessions from {acc2taxid_file}")
                log.info(f"Total unique accessions: {len(combined_df):,}")

                # Rebuild accession map from combined DataFrame
                self.acc_map = _build_accession_map_from_dataframe(combined_df, 'accession')

    @classmethod
    def load(cls, db_dir: str) -> 'TaxonomyDB':
        """
        Load pre-built taxonomy database from Parquet files.

        This is ~100x faster than building from dump files.

        Parameters
        ----------
        db_dir : str
            Directory containing Parquet files

        Returns
        -------
        TaxonomyDB
            Loaded taxonomy database
        """
        log.info(f"Loading taxonomy database from {db_dir}...")

        tax_db = TaxonomyDatabase.from_parquet(db_dir)

        # Load accession map if exists
        acc_map_file = os.path.join(db_dir, 'accession_map.parquet')
        acc_map = None
        if os.path.exists(acc_map_file):
            acc_map = AccessionMapping.from_parquet(acc_map_file)

        log.info("Database loaded successfully")
        return cls(tax_db, acc_map)

    def save(self, output_dir: str):
        """
        Save taxonomy database to Parquet files for instant loading.

        Parameters
        ----------
        output_dir : str
            Output directory path
        """
        log.info(f"Saving taxonomy database to {output_dir}...")

        os.makedirs(output_dir, exist_ok=True)

        # Save taxonomy tree
        self.tax_db.to_parquet(output_dir)

        # Save accession map if exists
        if self.acc_map is not None:
            acc_map_file = os.path.join(output_dir, 'accession_map.parquet')
            self.acc_map.to_parquet(acc_map_file)

        log.info("Database saved successfully")

    def get_name(self, taxid: int) -> Optional[str]:
        """
        Get scientific name for a taxid.

        Parameters
        ----------
        taxid : int
            Taxonomy ID

        Returns
        -------
        str or None
            Scientific name, or None if not found
        """
        return self.tax_db.get_name(taxid)

    def get_parent(self, taxid: int) -> Optional[int]:
        """
        Get parent taxid.

        Parameters
        ----------
        taxid : int
            Taxonomy ID

        Returns
        -------
        int or None
            Parent taxid, or None if not found
        """
        return self.tax_db.get_parent(taxid)

    def get_rank(self, taxid: int) -> Optional[str]:
        """
        Get rank name for a taxid.

        Parameters
        ----------
        taxid : int
            Taxonomy ID

        Returns
        -------
        str or None
            Rank name (e.g., 'species', 'genus'), or None if not found
        """
        return self.tax_db.get_rank(taxid)

    def get_lineage(self, taxid: int) -> Optional[List[int]]:
        """
        Get full lineage from root to taxid.

        Parameters
        ----------
        taxid : int
            Taxonomy ID

        Returns
        -------
        list of int or None
            List of taxids from root to the given taxid, or None if not found
        """
        return self.tax_db.get_lineage(taxid)

    def get_lineage_names(self, taxid: int) -> Optional[List[str]]:
        """
        Get full lineage names from root to taxid.

        Parameters
        ----------
        taxid : int
            Taxonomy ID

        Returns
        -------
        list of str or None
            List of names from root to the given taxid, or None if not found
        """
        lineage = self.get_lineage(taxid)
        if lineage is None:
            return None
        return [self.get_name(tid) for tid in lineage]

    def compute_lca(self, taxids: Union[int, List[int]], taxid2: Optional[int] = None) -> Optional[int]:
        """
        Compute Lowest Common Ancestor of taxids.

        Parameters
        ----------
        taxids : int or list of int
            Single taxid, pair of taxids, or list of taxids
        taxid2 : int, optional
            Second taxid if passing pair as two arguments

        Returns
        -------
        int or None
            LCA taxid, or None if not found

        Examples
        --------
        >>> tax.compute_lca(9606, 9598)  # Two taxids
        9604

        >>> tax.compute_lca([9606, 9598, 9593])  # List of taxids
        9604
        """
        if taxid2 is not None:
            # Called as compute_lca(taxid1, taxid2)
            return self.tax_db.compute_lca(taxids, taxid2)
        elif isinstance(taxids, (list, tuple)):
            # Called as compute_lca([taxid1, taxid2, ...])
            if len(taxids) == 0:
                return None
            elif len(taxids) == 1:
                return taxids[0]
            else:
                return self.tax_db.compute_lca_multi(list(taxids))
        else:
            # Single taxid
            return taxids

    def acc2taxid(self, accession: str) -> Optional[int]:
        """
        Get taxid for an accession.

        Parameters
        ----------
        accession : str
            Accession ID (e.g., 'NC_000001.11')

        Returns
        -------
        int or None
            Taxid if found, None otherwise
        """
        if self.acc_map is None:
            raise ValueError("Accession map not loaded")
        return self.acc_map.get_taxid(accession)

    def acc2taxid_batch(self, accessions: List[str]):
        """
        Get taxids for a batch of accessions (faster than individual queries).

        Parameters
        ----------
        accessions : list of str
            List of accession IDs

        Returns
        -------
        numpy.ndarray
            Array of taxids (int32), -1 for not found
        """
        if self.acc_map is None:
            raise ValueError("Accession map not loaded")
        return self.acc_map.get_taxids_batch(accessions)

    def build_lca_cache(self, taxids: List[int]):
        """
        Build LCA cache for O(1) queries on frequent taxid pairs.

        Parameters
        ----------
        taxids : list of int
            List of taxids to cache
        """
        self.tax_db.build_lca_cache(taxids)

    def get_taxid_from_reference(self, reference: str) -> Optional[int]:
        """
        Convenience method to get taxid from a reference name.

        This tries several common patterns:
        1. Direct lookup of reference
        2. Remove version suffix (e.g., '.1')
        3. Extract accession from concatenated reference (e.g., 'ACC1_ACC2' -> 'ACC1')

        Parameters
        ----------
        reference : str
            Reference name

        Returns
        -------
        int or None
            Taxid if found, None otherwise
        """
        if self.acc_map is None:
            raise ValueError("Accession map not loaded")

        # Try direct lookup
        taxid = self.acc2taxid(reference)
        if taxid is not None:
            return taxid

        # Try without version
        if '.' in reference:
            base = reference.rsplit('.', 1)[0]
            taxid = self.acc2taxid(base)
            if taxid is not None:
                return taxid

        # Try first part of concatenated reference
        if '_' in reference:
            first = reference.split('_')[0]
            taxid = self.acc2taxid(first)
            if taxid is not None:
                return taxid

        return None

    def get_taxonomy_dict(self, taxid: int) -> Optional[dict]:
        """
        Get full taxonomy information as a dictionary.

        Parameters
        ----------
        taxid : int
            Taxonomy ID

        Returns
        -------
        dict or None
            Dictionary with keys: taxid, name, rank, parent_taxid, lineage, lineage_names
        """
        name = self.get_name(taxid)
        if name is None:
            return None

        return {
            'taxid': taxid,
            'name': name,
            'rank': self.get_rank(taxid),
            'parent_taxid': self.get_parent(taxid),
            'lineage': self.get_lineage(taxid),
            'lineage_names': self.get_lineage_names(taxid),
        }

    @property
    def n_nodes(self) -> int:
        """Number of taxonomy nodes."""
        return self.tax_db.n_nodes

    @property
    def max_taxid(self) -> int:
        """Maximum taxid value."""
        return self.tax_db.max_taxid

    @property
    def n_accessions(self) -> int:
        """Number of accession mappings."""
        return self.acc_map.n_entries if self.acc_map is not None else 0

    def __repr__(self) -> str:
        return (
            f"TaxonomyDB(n_nodes={self.n_nodes:,}, "
            f"max_taxid={self.max_taxid:,}, "
            f"n_accessions={self.n_accessions:,})"
        )
