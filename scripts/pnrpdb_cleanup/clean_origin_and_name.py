#!/usr/bin/env python3
"""
clean_origin_and_name.py

- Default TSV path: ./pnrpdb2.tsv
- Simplify Origin by removing literal "from " occurrences (case-insensitive),
  then apply unicode normalization, mojibake fixes, dash/space normalization,
  and small exact replacements.
- Simplify Name according to duplication / plural-prefix rules.

Usage examples:
  python clean_origin_and_name.py
  python clean_origin_and_name.py --dry-run --show-origin-examples 20
  python clean_origin_and_name.py --path ./myfile.tsv --no-backup
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import re
import shutil
import sys
import unicodedata
from typing import List, Tuple

# -----------------------
# Configuration / Maps
# -----------------------

# Exact origin fixes (explicit one-off replacements)
ORIGIN_EXACT_FIXES = {
    "Streptomyces sp. strain G√∂-GS12": "Streptomyces sp. Gö-GS12",
    "Streptomyces sp. T√º 6071": "Streptomyces sp. Tü 6071",
    "Stigmatella aurantiaca DW4/3‚Äì1": "Stigmatella aurantiaca DW4/3-1",
    "Streptomyces sp. NRRL S-1868": "Streptomyces sp. NRRL S-1868",
    "Bacillus ¬† sp. SNA-60-367": "Bacillus sp. SNA-60-367",
    "Herbaspirillum seropedicae‚ÄÖZ67": "Herbaspirillum seropedicae Z67",
    "Lyngbya bouillonii PAL08‚Äì16": "Lyngbya bouillonii PAL08-16",
    "Streptomyces sp. (strain Tü 6071)": "Streptomyces sp. Tü 6071"
}

# Mojibake or recurring fragment fixes (extend this map as needed)
MOJIBAKE_FIXES = {
    "‚Äì": "-",   # common mojibake for en-dash
    "‚Äî": "-",   # another dash variant
    "‚Äö": " ",   # weird fragment -> space
    "‚Äö": " ",
    "‚ÄÖ": " ",
    "√∂": "ö",    # sometimes seen in mojibake
    "√º": "ü",
    "¬†": "",     # stray artifact -> remove
    "â€“": "-",   # another encoding artifact for dash
    "â€”": "-",   # em-dash variant
}

# Regexes
FROM_LITERAL_RE = re.compile(r"\bfrom\s+", flags=re.IGNORECASE)  # remove literal "from " occurrences
DASHES_RE = re.compile(r"[‐\u2010-\u2015\u2212-]")  # many dash-like characters
ODD_SPACES_RE = re.compile(r"[\u00A0\u2007\u202F\u2009]")  # non-breaking and narrow spaces

# -----------------------
# Origin simplification
# -----------------------


def remove_from_literal(s: str) -> str:
    """Remove all occurrences of literal 'from ' (case-insensitive)."""
    return FROM_LITERAL_RE.sub("", s)


def apply_mojibake_and_norm(s: str) -> str:
    """
    Apply exact fixes, mojibake replacements and normalize unicode/dashes/spaces.

    Steps:
      1. If exact-match exists in ORIGIN_EXACT_FIXES -> return that immediately.
      2. Unicode normalize (NFC)
      3. Replace known mojibake fragments.
      4. Normalize dash characters to ASCII '-'
      5. Replace odd spaces with normal space, collapse multiple spaces, strip.
    """
    if s in ORIGIN_EXACT_FIXES:
        return ORIGIN_EXACT_FIXES[s]

    # Unicode normalization
    s2 = unicodedata.normalize("NFC", s)

    # Mojibake fragment fixes
    for bad, good in MOJIBAKE_FIXES.items():
        if bad in s2:
            s2 = s2.replace(bad, good)

    # Normalize dashes
    s2 = DASHES_RE.sub("-", s2)

    # Normalize odd spaces to normal space
    s2 = ODD_SPACES_RE.sub(" ", s2)

    # Collapse multiple spaces, strip
    s2 = re.sub(r"\s+", " ", s2).strip()

    return s2


def simplify_origin(origin: str) -> Tuple[str, bool]:
    """
    Full origin simplification:
      - Remove literal "from " occurrences
      - Then apply mojibake fixes & normalization
    Returns (new_origin, changed_bool)
    """
    if origin is None:
        return origin, False

    s = str(origin)
    if not s:
        return origin, False

    # 1) remove "from " occurrences
    removed = remove_from_literal(s)

    # 2) apply normalization + mojibake fixes
    cleaned = apply_mojibake_and_norm(removed)

    changed = cleaned != s
    return cleaned, changed

#### OLD VERSION (more tricky "From" removal)
"""

FROM_RE = re.compile(r"\bfrom\b", flags=re.IGNORECASE)
LEADING_ARTICLE_RE = re.compile(r"^(?:a|the)\s+", flags=re.IGNORECASE)


def simplify_origin(origin: str) -> Tuple[str, bool]:
    
    # Returns (new_origin, changed).
    
    if origin is None:
        return origin, False

    s = str(origin).strip()
    if not s:
        return origin, False

    m = FROM_RE.search(s)
    if not m:
        return origin, False

    # Keep everything after the first "from"
    tail = s[m.end():].strip()

    # Drop leading "a " / "the "
    tail2 = LEADING_ARTICLE_RE.sub("", tail).strip()

    # If tail becomes empty, keep original (defensive)
    if not tail2:
        return origin, False

    return tail2, (tail2 != origin)


"""


# --- PNPdb-specific name cleanup ---

PNPDB_ID_RE = re.compile(r"^PNPdb", flags=0)

E_QUOTED_NAME_RE = re.compile(r"^E'(.+)'$")

def normalize_pnpdb_name(name: str) -> str:
    """
    Normalize Name strings for PNPdb entries (ID starts with 'PNPdb').

    - Remove all double-quote characters (") anywhere in the string.
      Examples:
        \"""SB_245255;Cyclo...\""" -> SB_245255;Cyclo...
        \"""Pyoverdin_\"""\"""1b\"""\"""\"""
          -> Pyoverdin_1b

    - Convert templated E'name' -> name
      Examples:
        E'actinoramide_B' -> actinoramide_B
        E'marformycin_E'  -> marformycin_E

    Trims surrounding whitespace after transformations.
    """
    if name is None:
        return ""

    s = str(name).strip()

    # E'name' template
    m = E_QUOTED_NAME_RE.match(s)
    if m:
        s = m.group(1)

    # Remove all double quotes (Numbers/CSV export artifacts, etc.)
    if '"' in s:
        s = s.replace('"', "")

    return s.strip()

# -----------------------
# Name simplification
# -----------------------


def simplify_name(name: str) -> Tuple[str, bool]:
    """
    Simplify duplicated / plural-prefixed names.

    Rules:
      - Split on '_'
      - If there exists a split where right starts with left:
          * if right == left -> keep left
          * else -> keep right (the more specific second half)
      - Otherwise, if first part endswith 's' and first[:-1] == second:
          -> drop first (keep rest)
    """
    if name is None:
        return name, False

    s = str(name).strip()
    if not s:
        return name, False

    parts = s.split("_")
    if len(parts) < 2:
        return name, False

    new_parts = None

    # Try duplication / extension pattern
    for i in range(1, len(parts)):
        left = parts[:i]
        right = parts[i:]
        if len(right) >= len(left) and right[: len(left)] == left:
            # found a match
            if right == left:
                new_parts = left  # exact duplication
            else:
                new_parts = right  # extension, keep more specific
            break

    # Try plural-prefix rule
    if new_parts is None:
        first = parts[0]
        if first.endswith("s") and len(parts) >= 2:
            singular = first[:-1]
            if singular == parts[1]:
                new_parts = parts[1:]

    if new_parts is None:
        return name, False

    new_name = "_".join(new_parts)
    if new_name == s:
        return name, False

    return new_name, True


# -----------------------
# TSV utilities
# -----------------------


def read_tsv(path: str) -> Tuple[List[str], List[dict]]:
    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("TSV appears to have no header row.")
        rows = list(reader)
        return reader.fieldnames, rows


def write_tsv(path: str, fieldnames: List[str], rows: List[dict]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            quoting=csv.QUOTE_MINIMAL,
        )
        writer.writeheader()
        for r in rows:
            out = {k: r.get(k, "") for k in fieldnames}
            writer.writerow(out)


def backup_file(path: str) -> str:
    ts = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = f"{path}.bak_{ts}"
    shutil.copy2(path, backup_path)
    return backup_path


# -----------------------
# Main
# -----------------------


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--path",
        default="./pnrpdb2.tsv",
        help="Path to the TSV to modify in place (default: ./pnrpdb2.tsv).",
    )

    ap.add_argument(
        "--id-col",
        default="ID",
        help='Column name for ID (used for PNPdb-specific Name cleanup) (default: "ID").',
    )
    ap.add_argument(
        "--origin-col",
        default="Origin",
        help='Column name to edit for origin simplification (default: "Origin").',
    )
    ap.add_argument(
        "--name-col",
        default="Name",
        help='Column name to edit for name simplification (default: "Name").',
    )
    ap.add_argument(
        "--no-backup",
        action="store_true",
        help="Do not create a timestamped .bak_* backup before writing.",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Do not write anything; just report how many rows would change.",
    )
    ap.add_argument(
        "--show-origin-examples",
        type=int,
        default=10,
        help="Print up to N before/after examples of changed Origin rows (default: 10).",
    )
    ap.add_argument(
        "--show-name-examples",
        type=int,
        default=10,
        help="Print up to N before/after examples of changed Name rows (default: 10).",
    )
    args = ap.parse_args()

    path = args.path
    if not os.path.exists(path):
        print(f"ERROR: file not found: {path}", file=sys.stderr)
        return 2

    fieldnames, rows = read_tsv(path)

    if args.origin_col not in fieldnames:
        print(
            f'WARNING: origin column "{args.origin_col}" not found. Origin will not be modified.',
            file=sys.stderr,
        )

    if args.name_col not in fieldnames:
        print(
            f'WARNING: name column "{args.name_col}" not found. Name will not be modified.',
            file=sys.stderr,
        )

    if args.id_col not in fieldnames:
        print(
            f'WARNING: id column "{args.id_col}" not found. PNPdb-specific Name cleanup will be skipped.',
            file=sys.stderr,
        )

    origin_changed = 0
    name_changed = 0
    origin_examples: List[Tuple[str, str]] = []
    name_examples: List[Tuple[str, str]] = []

    for row in rows:
        # Origin
        if args.origin_col in fieldnames:
            old_o = row.get(args.origin_col, "")
            new_o, changed_o = simplify_origin(old_o)
            if changed_o:
                row[args.origin_col] = new_o
                origin_changed += 1
                if len(origin_examples) < args.show_origin_examples:
                    origin_examples.append((old_o, new_o))

        # Name
        if args.name_col in fieldnames:
            orig_n = row.get(args.name_col, "") or ""
            proc_n = orig_n  # will hold value after normalization + simplification

            # PNPdb-specific cleanup: remove stray quotes and E'...' template
            id_val = str(row.get(args.id_col, "") or "")
            normalized_changed = False
            if PNPDB_ID_RE.match(id_val):
                normalized = normalize_pnpdb_name(orig_n)
                if normalized != orig_n:
                    normalized_changed = True
                    proc_n = normalized
            # If not PNPdb or nothing changed, proc_n is still orig_n

            # Now apply general simplification (dedup / plural rules)
            simplified_n, simplified_changed = simplify_name(proc_n)

            # Determine final value and whether anything changed at all
            final_n = simplified_n if simplified_changed else proc_n
            any_changed = normalized_changed or simplified_changed or (final_n != orig_n)

            if any_changed:
                row[args.name_col] = final_n
                name_changed += 1
                # record example using original (before any modification) -> final
                if len(name_examples) < args.show_name_examples:
                    name_examples.append((orig_n, final_n))

    print(f"Rows total: {len(rows)}")
    if args.origin_col in fieldnames:
        print(f'Rows with changed "{args.origin_col}": {origin_changed}')
    if args.name_col in fieldnames:
        print(f'Rows with changed "{args.name_col}": {name_changed}')

    if origin_examples:
        print("\nOrigin examples (before -> after):")
        for b, a in origin_examples:
            print(f"- {b!r} -> {a!r}")

    if name_examples:
        print("\nName examples (before -> after):")
        for b, a in name_examples:
            print(f"- {b!r} -> {a!r}")

    if args.dry_run:
        print("\nDry-run: no files written.")
        return 0

    if not args.no_backup:
        bpath = backup_file(path)
        print(f"\nBackup written: {bpath}")

    write_tsv(path, fieldnames, rows)
    print(f"Updated TSV written in place: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())