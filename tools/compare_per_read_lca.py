#!/usr/bin/env python3
import argparse
import gzip
import re
from typing import Dict, Tuple

# metaDMG: columns: queryid, seq, len, naln, gc, lca, taxa_path
#   lca looks like: 79494:"s__SOJZ01 sp004376295":"species"
MD_LCA_RE = re.compile(r"^(\d+):\".*\":\"([a-z_]+)\"$")


def load_metadmg(path: str) -> Dict[str, Tuple[int, str, int]]:
    m: Dict[str, Tuple[int, str, int]] = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="ignore") as fh:
        header_skipped = False
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            if not header_skipped:
                # expect a header line with columns
                if line.lower().startswith("queryid"):
                    header_skipped = True
                    continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            rid = parts[0]
            try:
                naln = int(parts[3])
            except ValueError:
                naln = -1
            lca_field = parts[5]
            mobj = MD_LCA_RE.match(lca_field)
            if not mobj:
                # sometimes quotes or escapes can differ; try to salvage taxid as leading integer
                taxid_str = lca_field.split(":", 1)[0]
                try:
                    taxid = int(taxid_str)
                except Exception:
                    continue
                rank = ""
            else:
                taxid = int(mobj.group(1))
                rank = mobj.group(2)
            m[rid] = (taxid, rank, naln)
    return m


def load_ours(path: str) -> Dict[str, Tuple[int, str, int]]:
    o: Dict[str, Tuple[int, str, int]] = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="ignore") as fh:
        # header
        fh.readline()
        # expected: read_name\tlca_taxid\tlca_rank\tn_aln\ttax_path\tnorm_ref\tnorm_ref_len
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            rid = parts[0]
            try:
                taxid = int(parts[1])
            except ValueError:
                continue
            rank = parts[2]
            try:
                naln = int(parts[3])
            except ValueError:
                naln = -1
            o[rid] = (taxid, rank, naln)
    return o


def main():
    ap = argparse.ArgumentParser(
        description="Compare per-read LCA outputs (metaDMG vs our per-read)"
    )
    ap.add_argument(
        "--metadmg",
        required=True,
        help="Path to metaDMG per-read LCA (e.g., test-md.lca.gz)",
    )
    ap.add_argument(
        "--ours",
        required=True,
        help="Path to our per-read TSV (e.g., test.strict.per-read.v3.tsv.gz)",
    )
    ap.add_argument(
        "--show", type=int, default=10, help="Show first N mismatches (default 10)"
    )
    args = ap.parse_args()

    md = load_metadmg(args.metadmg)
    print(f"metaDMG reads: {len(md):,}")
    ours = load_ours(args.ours)
    print(f"Our reads:    {len(ours):,}")

    common = set(md.keys()) & set(ours.keys())
    only_md = set(md.keys()) - set(ours.keys())
    only_ours = set(ours.keys()) - set(md.keys())

    match = 0
    mismatches = []
    rank_diff = 0
    for rid in common:
        t_md, r_md, _ = md[rid]
        t_ou, r_ou, _ = ours[rid]
        if t_md == t_ou:
            match += 1
        else:
            mismatches.append((rid, t_md, r_md, t_ou, r_ou))
        if r_md != r_ou:
            rank_diff += 1

    print(f"Common reads: {len(common):,}")
    print(f"Matches:      {match:,}")
    print(f"Mismatches:   {len(mismatches):,}")
    print(f"Only metaDMG: {len(only_md):,}")
    print(f"Only ours:    {len(only_ours):,}")
    print(f"Rank diffs among common: {rank_diff:,}")

    if mismatches:
        print("\nFirst mismatches:")
        for row in mismatches[: args.show]:
            rid, t_md, r_md, t_ou, r_ou = row
            print(f"{rid}\tmd:{t_md}({r_md})\tours:{t_ou}({r_ou})")


if __name__ == "__main__":
    main()
