#!/usr/bin/env python3
"""
Quick checker for accession -> taxid mapping inside a Parquet taxonomy DB,
and optional cross-check against a gzipped acc2taxid map (e.g., metaDMG input).

Usage examples:
  - Check a single accession against a Parquet taxonomy DB:
      python tools/check_accession_mapping.py --taxonomy-db /path/to/db --accession GCF_001514235.1

  - Check a list of accessions from a file and compare to acc2tax gz:
      python tools/check_accession_mapping.py \
        --taxonomy-db /path/to/db \
        --accessions-file acc_list.txt \
        --acc2tax-gz /path/to/acc2taxid.map.gz

This script uses DuckDB to read Parquet files (no pyarrow dependency required).
"""
import argparse
import gzip
import os
import sys
from typing import Dict, Iterable, List, Optional, Tuple

import duckdb


def _norm_candidates(acc: str) -> List[str]:
    """Return candidate accession strings to try (exact then version-stripped)."""
    cands = [acc]
    if acc.startswith(("GCF_", "GCA_")) and "." in acc:
        base = acc.split(".", 1)[0]
        if base and base != acc:
            cands.append(base)
    return cands


def _query_parquet_taxid(
    conn: duckdb.DuckDBPyConnection, parquet_path: str, accession: str
) -> Optional[int]:
    """Query accession_map.parquet for a taxid by accession with simple fallbacks.

    Returns the first taxid found or None.
    """
    for cand in _norm_candidates(accession):
        try:
            res = conn.execute(
                """
                SELECT taxid
                FROM read_parquet(?)
                WHERE accession = ?
                LIMIT 1
                """,
                [parquet_path, cand],
            ).fetchone()
        except Exception as e:
            print(f"[err] DuckDB query failed for {cand}: {e}", file=sys.stderr)
            continue
        if res is not None and len(res) > 0 and res[0] is not None:
            return int(res[0])
    return None


def _load_node(
    conn: duckdb.DuckDBPyConnection, nodes_parquet: str, taxid: int
) -> Optional[Tuple[int, int, str, str]]:
    """Return (taxid, parent_taxid, rank, name) for a given taxid or None."""
    try:
        row = conn.execute(
            """
            SELECT taxid, parent_taxid, rank, name
            FROM read_parquet(?)
            WHERE taxid = ?
            LIMIT 1
            """,
            [nodes_parquet, taxid],
        ).fetchone()
    except Exception as e:
        print(f"[err] DuckDB query failed for taxid {taxid}: {e}", file=sys.stderr)
        return None
    if row:
        return int(row[0]), int(row[1]), str(row[2]), str(row[3])
    return None


def _build_lineage(
    conn: duckdb.DuckDBPyConnection, nodes_parquet: str, taxid: int
) -> List[Tuple[int, str, str]]:
    """Ascend parent_taxid until root; returns list of (taxid, rank, name) from root->node."""
    lineage: List[Tuple[int, str, str]] = []
    current = taxid
    seen = set()
    while True:
        node = _load_node(conn, nodes_parquet, current)
        if not node:
            break
        tid, parent_tid, rank, name = node
        lineage.append((tid, rank, name))
        if tid == parent_tid:
            break  # root
        if tid in seen or parent_tid <= 0:
            break
        seen.add(tid)
        current = parent_tid
    lineage.reverse()
    return lineage


def _read_acc2tax_gz(acc2tax_gz: str, targets: Iterable[str]) -> Dict[str, int]:
    """Return mapping accession->taxid from a gz TSV for target accessions only.

    The gz file is expected to have at least two columns: accession and taxid.
    We'll scan linearly and stop after all targets are found.
    """
    targets_set = set(targets)
    found: Dict[str, int] = {}
    with gzip.open(acc2tax_gz, "rt", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not targets_set:
                break
            if not line or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            acc, tax = parts[0], parts[1]
            if acc in targets_set:
                try:
                    found[acc] = int(tax)
                except ValueError:
                    continue
                targets_set.remove(acc)
    return found


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Check accession->taxid mapping in Parquet taxonomy DB"
    )
    ap.add_argument(
        "--taxonomy-db",
        required=True,
        help="Directory with nodes.parquet and accession_map.parquet",
    )
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument(
        "--accession",
        action="append",
        help="Accession to check (can be given multiple times)",
    )
    grp.add_argument("--accessions-file", help="File with one accession per line")
    ap.add_argument(
        "--acc2tax-gz", help="Optional path to gzipped acc2taxid map to compare against"
    )
    args = ap.parse_args()

    db_dir = args.taxonomy_db
    nodes_parquet = os.path.join(db_dir, "nodes.parquet")
    accmap_parquet = os.path.join(db_dir, "accession_map.parquet")

    if not os.path.isdir(db_dir):
        print(f"[fatal] taxonomy-db directory not found: {db_dir}", file=sys.stderr)
        return 2
    for p in (nodes_parquet, accmap_parquet):
        if not os.path.exists(p):
            print(f"[fatal] required Parquet file missing: {p}", file=sys.stderr)
            return 2

    # Build accession list
    accessions: List[str] = []
    if args.accession:
        for a in args.accession:
            accessions.append(a.strip())
    if args.accessions_file:
        with open(args.accessions_file, "r", encoding="utf-8") as fh:
            for line in fh:
                a = line.strip()
                if a:
                    accessions.append(a)

    if not accessions:
        print("[fatal] no accessions provided", file=sys.stderr)
        return 2

    # Optional acc2tax gz check
    acc2tax_map: Dict[str, int] = {}
    if args.acc2tax_gz:
        if not os.path.exists(args.acc2tax_gz):
            print(f"[warn] acc2tax gz not found: {args.acc2tax_gz}", file=sys.stderr)
        else:
            acc2tax_map = _read_acc2tax_gz(args.acc2tax_gz, accessions)

    conn = duckdb.connect()
    print(f"Using taxonomy DB: {db_dir}")
    print()

    for acc in accessions:
        print(f"=== {acc} ===")
        parquet_taxid = _query_parquet_taxid(conn, accmap_parquet, acc)
        if parquet_taxid is None:
            print("Parquet: NOT FOUND")
        else:
            node = _load_node(conn, nodes_parquet, parquet_taxid)
            if node:
                tid, parent_tid, rank, name = node
                lineage = _build_lineage(conn, nodes_parquet, tid)
                lineage_str = ";".join([f"{r}__{n}" for (_, r, n) in lineage])
                print(f"Parquet: taxid={tid} rank={rank} name={name}")
                print(f"         lineage={lineage_str}")
            else:
                print(
                    f"Parquet: taxid={parquet_taxid} (node not found in nodes.parquet)"
                )

        if acc in acc2tax_map:
            print(f"metaDMG acc2tax: taxid={acc2tax_map[acc]}")
        elif args.acc2tax_gz:
            print("metaDMG acc2tax: NOT FOUND (in provided gz)")

        print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
