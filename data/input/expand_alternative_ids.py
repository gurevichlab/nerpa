#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Iterable, Sequence


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            'Expand a table by adding one row per Alternative ID. '
            'Each expanded row has ID=alt_id and Alternative IDs set to the other IDs in the group.'
        )
    )
    p.add_argument('input', type=Path, help='Input file (TSV by default).')
    p.add_argument('-o', '--output', type=Path, default=None, help='Output path (default: stdout).')
    p.add_argument(
        '--delimiter',
        default='\t',
        help=r"CSV delimiter (default: '\t').",
    )
    p.add_argument('--id-col', default='ID', help="ID column name (default: 'ID').")
    p.add_argument(
        '--alt-col',
        default='Alternative IDs',
        help="Alternative IDs column name (default: 'Alternative IDs').",
    )
    p.add_argument(
        '--alt-sep',
        default=';',
        help="Separator to use when writing Alternative IDs in expanded rows (default: ';').",
    )
    return p.parse_args(argv)


def _split_alt_ids(s: str) -> list[str]:
    # Accept either ';' or ',' in the input; trim whitespace; drop empties.
    parts: list[str] = []
    for chunk in s.replace(',', ';').split(';'):
        item = chunk.strip()
        if item:
            parts.append(item)
    return parts


def iter_expanded_rows(
    rows: Iterable[dict[str, str]],
    *,
    id_col: str,
    alt_col: str,
    alt_sep: str,
) -> Iterable[dict[str, str]]:
    for row in rows:
        yield row  # keep original

        orig_id: str = row.get(id_col, '')
        alt_raw: str = row.get(alt_col, '')
        alt_ids: list[str] = _split_alt_ids(alt_raw)

        if not alt_ids:
            continue

        group_ids: list[str] = [orig_id, *alt_ids]

        for alt_id in alt_ids:
            new_row: dict[str, str] = dict(row)
            new_row[id_col] = alt_id

            others: list[str] = [x for x in group_ids if x and x != alt_id]
            new_row[alt_col] = alt_sep.join(others)

            yield new_row


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)

    inp = args.input.open('r', newline='', encoding='utf-8')
    try:
        reader = csv.DictReader(inp, delimiter=args.delimiter)
        if reader.fieldnames is None:
            raise SystemExit('Input has no header row.')

        out_fh = sys.stdout if args.output is None else args.output.open('w', newline='', encoding='utf-8')
        try:
            writer = csv.DictWriter(out_fh, fieldnames=reader.fieldnames, delimiter=args.delimiter)
            writer.writeheader()

            for out_row in iter_expanded_rows(
                reader,
                id_col=args.id_col,
                alt_col=args.alt_col,
                alt_sep=args.alt_sep,
            ):
                writer.writerow(out_row)
        finally:
            if out_fh is not sys.stdout:
                out_fh.close()
    finally:
        inp.close()

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
