#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from src.build_output.chem_helper import MolRecord
from src.rban_parsing.rban_parser import Parsed_rBAN_Record


def build_id_to_smiles(records_by_id: dict[str, Any]) -> dict[str, str]:
    id_to_smiles: dict[str, str] = {}

    for rec_id, rec_dict in records_by_id.items():
        record: Parsed_rBAN_Record = Parsed_rBAN_Record.from_dict(rec_dict)
        smiles: str = MolRecord.from_rban_record(record).to_canonical_smiles()
        id_to_smiles[str(rec_id)] = smiles

    return id_to_smiles


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            'Read a JSON dict of Parsed_rBAN_Record objects keyed by id and output a JSON dict '
            'mapping id -> canonical SMILES.'
        )
    )
    p.add_argument(
        'input',
        type=Path,
        help='Path to input JSON file (a dict keyed by record id)',
    )
    p.add_argument(
        '-o',
        '--output',
        type=Path,
        help='Write output JSON here',
    )
    p.add_argument(
        '--indent',
        type=int,
        default=2,
        help='Pretty-print JSON with this indent level (default: 2). Use 0 for compact.',
    )
    return p.parse_args()


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding='utf-8'))


def write_json(data: Any, path: Path | None, indent: int) -> None:
    import sys

    out_str: str
    if indent <= 0:
        out_str = json.dumps(data, ensure_ascii=False, separators=(',', ':'))
    else:
        out_str = json.dumps(data, ensure_ascii=False, indent=indent)

    if path is None:
        sys.stdout.write(out_str)
        sys.stdout.write('\n')
        return

    path.write_text(out_str + '\n', encoding='utf-8')


def main() -> None:
    args = parse_args()

    raw: Any = read_json(args.input_json)
    if not isinstance(raw, dict):
        raise SystemExit('Input JSON must be a dict keyed by record id.')

    id_to_smiles = build_id_to_smiles(raw)
    write_json(id_to_smiles, args.output, args.indent)


if __name__ == '__main__':
    main()
