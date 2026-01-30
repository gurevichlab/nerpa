#!/usr/bin/env python3
"""
Collapse pnrpdb2.tsv by unique InChiKey.

Input:  ./pnrpdb2.tsv
Output: ./pnrpdb2_unique_InChiKey.tsv

Rules (summary):
- One row per InChiKey.
- Select primary row by ID prefix priority (BGC00 > NOR > NPA > PNP), then lexicographically.
- Name: dedupe case-/separator-insensitively, preserve first spelling; drop names
        whose key is a strict prefix of another (keep longer).
- SMILES: from primary
- Origin: unique list starting from primary, joined by "; "
- Source: unique list starting from primary, joined by "; "
- InChiKey: from primary
- Alternative IDs: list of other IDs in priority order joined by "; "
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import unicodedata
import sys
from collections import defaultdict
from typing import Any, Dict, List, Tuple

# Priority order of ID prefixes (lower index = higher priority)
ID_PREFIX_PRIORITY = ["BGC00", "NOR", "NPA", "PNP"]

ALTERNATIVE_IDS_COLNAME = "Alternative IDs"


def id_priority(id_value: str) -> Tuple[int, str]:
    """
    Return a (priority_rank, id_value) tuple for sorting.

    - priority_rank: index into ID_PREFIX_PRIORITY if matched,
                     len(ID_PREFIX_PRIORITY) if not matched (lowest priority).
    - id_value: used as tie-breaker lexicographically.
    """
    if id_value is None:
        return len(ID_PREFIX_PRIORITY), ""

    s = str(id_value)
    for rank, prefix in enumerate(ID_PREFIX_PRIORITY):
        if s.startswith(prefix):
            return rank, s

    return len(ID_PREFIX_PRIORITY), s


def unique_preserve_order(values: List[str]) -> List[str]:
    """Return list of unique strings, preserving first occurrence order."""
    seen = set()
    out = []
    for v in values:
        if v is None:
            continue
        v = str(v)
        if v == "":
            continue
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out


# -----------------------
# Name canonicalization helpers
# -----------------------

# treat various separators (space, underscore, hyphen, many unicode dashes) as equivalent
_SEP_RE = re.compile(r"[\s_\-‐-‒–—−]+")
# remove any non-alphanumeric after folding
_NONALNUM_RE = re.compile(r"[^0-9a-z]+")


def _name_compare_key(name: str) -> str:
    """
    Comparison key for names:
    - Unicode normalize NFC
    - casefold (better than lower)
    - make 'α' equivalent to 'alpha'
    - remove separators (space, underscore, hyphen, unicode dashes)
    - remove remaining non-alphanumeric chars
    """
    if not name:
        return ""
    s = unicodedata.normalize("NFC", name).strip()
    s = s.casefold()
    # map Greek alpha to literal "alpha"
    s = s.replace("α", "alpha")
    # remove separators
    s = _SEP_RE.sub("", s)
    # drop remaining non-alphanumeric
    s = _NONALNUM_RE.sub("", s)
    return s


def collapse_names_priority_prefix(raw_names_in_priority_order: List[str]) -> str:
    """
    Deduplicate names in priority order using comparison keys (casefold+separator-insensitive).
    Then drop any name whose key is a strict prefix of another key (keep the longer).
    Preserve original spellings for the kept names, in the order of priority.
    Return joined string with "; ".
    """
    # clean empties
    raw = [n for n in (raw_names_in_priority_order or []) if (n or "").strip() != ""]
    if not raw:
        return ""

    # Step 1: dedupe by comparison key, keep first spelling encountered
    kept_names: List[str] = []
    kept_keys: List[str] = []
    seen_keys = set()
    for n in raw:
        k = _name_compare_key(n)
        if not k:
            continue
        if k in seen_keys:
            continue
        seen_keys.add(k)
        kept_names.append(n)
        kept_keys.append(k)

    if len(kept_names) <= 1:
        return kept_names[0] if kept_names else ""

    # Step 2: prefix pruning: if key_i is strict prefix of key_j, drop i
    drop = [False] * len(kept_names)
    for i, ki in enumerate(kept_keys):
        if not ki:
            continue
        for j, kj in enumerate(kept_keys):
            if i == j:
                continue
            if kj and ki != kj and kj.startswith(ki):
                drop[i] = True
                break

    final = [n for n, d in zip(kept_names, drop) if not d]
    return "; ".join(final)


# -----------------------
# Collapse logic
# -----------------------


def collapse_group(
    rows: List[Dict[str, Any]],
    id_col: str,
    name_col: str,
    smiles_col: str,
    origin_col: str,
    inchikey_col: str,
    source_col: str,
) -> Dict[str, Any]:
    """
    Collapse a list of rows (all having the same InChiKey) into a single row.
    """
    if not rows:
        raise ValueError("collapse_group called with empty rows list")

    # Sort rows by ID priority (primary first)
    sorted_rows = sorted(rows, key=lambda r: id_priority(r.get(id_col, "")))

    primary = sorted_rows[0]
    collapsed = dict(primary)  # start from primary row

    # --- Name ---
    raw_names = [str(primary.get(name_col, "") or "")]
    for r in sorted_rows[1:]:
        raw_names.append(str(r.get(name_col, "") or ""))

    name_merged = collapse_names_priority_prefix(raw_names)
    collapsed[name_col] = name_merged or primary.get(name_col, "")

    # --- SMILES ---
    collapsed[smiles_col] = primary.get(smiles_col, "")

    # --- Origin ---
    origin_values = [primary.get(origin_col, "")]
    for r in sorted_rows[1:]:
        origin_values.append(r.get(origin_col, ""))
    origin_uniq = unique_preserve_order(origin_values)
    collapsed[origin_col] = "; ".join(origin_uniq) if origin_uniq else primary.get(origin_col, "")

    # --- InChiKey ---
    collapsed[inchikey_col] = primary.get(inchikey_col, "")

    # --- Source ---
    source_values = [primary.get(source_col, "")]
    for r in sorted_rows[1:]:
        source_values.append(r.get(source_col, ""))
    source_uniq = unique_preserve_order(source_values)
    collapsed[source_col] = "; ".join(source_uniq) if source_uniq else primary.get(source_col, "")

    # --- Alternative IDs ---
    alt_ids = []
    for r in sorted_rows[1:]:
        alt_id = r.get(id_col, "")
        if alt_id:
            alt_ids.append(str(alt_id))
    alt_ids_uniq = unique_preserve_order(alt_ids)
    collapsed[ALTERNATIVE_IDS_COLNAME] = "; ".join(alt_ids_uniq) if alt_ids_uniq else ""

    return collapsed


def main() -> int:
    parser = argparse.ArgumentParser(description="Collapse pnrpdb2.tsv by InChiKey into unique entries.")
    parser.add_argument("--input", default="./pnrpdb2.tsv", help='Input TSV path (default: "./pnrpdb2.tsv").')
    parser.add_argument(
        "--output", default="./pnrpdb2_unique_InChiKey.tsv", help='Output TSV path (default: "./pnrpdb2_unique_InChiKey.tsv").'
    )
    parser.add_argument("--id-col", default="ID", help='Column name for ID (default: "ID").')
    parser.add_argument("--name-col", default="Name", help='Column name for Name (default: "Name").')
    parser.add_argument("--smiles-col", default="SMILES", help='Column name for SMILES (default: "SMILES").')
    parser.add_argument("--origin-col", default="Origin", help='Column name for Origin (default: "Origin").')
    parser.add_argument("--inchikey-col", default="InChiKey", help='Column name for InChiKey (default: "InChiKey").')
    parser.add_argument("--source-col", default="Source", help='Column name for Source (default: "Source").')
    args = parser.parse_args()

    in_path = args.input
    out_path = args.output

    if not os.path.exists(in_path):
        print(f"ERROR: input file not found: {in_path}", file=sys.stderr)
        return 1

    # Read input TSV
    with open(in_path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            print("ERROR: input TSV has no header row.", file=sys.stderr)
            return 1
        fieldnames = list(reader.fieldnames)
        rows = list(reader)

    # Check required columns
    required_cols = [args.id_col, args.name_col, args.smiles_col, args.origin_col, args.inchikey_col, args.source_col]
    missing = [c for c in required_cols if c not in fieldnames]
    if missing:
        print("ERROR: missing required columns in TSV: " + ", ".join(missing), file=sys.stderr)
        return 1

    # Add Alternative IDs column if not present
    if ALTERNATIVE_IDS_COLNAME not in fieldnames:
        fieldnames.append(ALTERNATIVE_IDS_COLNAME)

    # Group by InChiKey
    groups: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    inchikey_col = args.inchikey_col
    for r in rows:
        key = r.get(inchikey_col, "")
        groups[key].append(r)

    # Collapse each group
    collapsed_rows: List[Dict[str, Any]] = []
    for inchikey, group_rows in groups.items():
        if len(group_rows) == 1:
            r = dict(group_rows[0])
            r.setdefault(ALTERNATIVE_IDS_COLNAME, "")
            collapsed_rows.append(r)
        else:
            collapsed = collapse_group(
                group_rows,
                id_col=args.id_col,
                name_col=args.name_col,
                smiles_col=args.smiles_col,
                origin_col=args.origin_col,
                inchikey_col=args.inchikey_col,
                source_col=args.source_col,
            )
            collapsed_rows.append(collapsed)

    # Sort for reproducibility (by selected ID then InChiKey)
    collapsed_rows.sort(key=lambda r: (str(r.get(args.id_col, "")), str(r.get(args.inchikey_col, ""))))

    # Write output TSV
    with open(out_path, "w", encoding="utf-8", newline="") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        for r in collapsed_rows:
            row_out = {k: r.get(k, "") for k in fieldnames}
            writer.writerow(row_out)

    print(f"Done. Collapsed {len(rows)} input rows into {len(collapsed_rows)} rows by InChiKey.")
    print(f"Output written to {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())