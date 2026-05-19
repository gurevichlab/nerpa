#!/usr/bin/env python3
from __future__ import annotations

import argparse
import yaml
from pathlib import Path
from typing import Any, List

from src.build_output.chem_helper import MolRecord
from src.rban_parsing.rban_parser import Parsed_rBAN_Record


def build_id_to_smiles(records_raw: List[Any]) -> dict[str, str]:
    id_to_smiles: dict[str, str] = {}

    for rec_dict in records_raw:
        rec_id: str = rec_dict['compound_id']
        record: Parsed_rBAN_Record = Parsed_rBAN_Record.from_dict(rec_dict)
        smiles: str = MolRecord.from_rban_record(record).to_canonical_smiles()
        id_to_smiles[str(rec_id)] = smiles

    return id_to_smiles


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            'Read a YAML list of Parsed_rBAN_Record objects and output a YAML dict '
            'mapping id -> canonical SMILES.'
        )
    )
    p.add_argument(
        'input',
        type=Path,
        help='Path to input YAML file (a dict keyed by record id)',
    )
    p.add_argument(
        '-o',
        '--output',
        type=Path,
        help='Write output YAML here',
    )
    p.add_argument(
        '--indent',
        type=int,
        default=2,
        help='Pretty-print YAML with this indent level (default: 2). Use 0 for compact.',
    )
    return p.parse_args()


def read_yaml(path: Path) -> Any:
    with path.open(encoding='utf-8') as f:
        return yaml.safe_load(f)


def write_yaml(data: Any, path: Path | None, indent: int) -> None:
    import sys

    out_str: str
    if indent <= 0:
        out_str = yaml.dump(data, default_flow_style=False, allow_unicode=True)
    else:
        out_str = yaml.dump(data, default_flow_style=False, allow_unicode=True, indent=indent)

    if path is None:
        sys.stdout.write(out_str)
        sys.stdout.write('\n')
        return

    path.write_text(out_str + '\n', encoding='utf-8')


def main() -> None:
    args = parse_args()

    raw: Any = read_yaml(args.input)
    if not isinstance(raw, list):
        raise SystemExit('Input YAML must be a list.')

    id_to_smiles = build_id_to_smiles(raw)
    write_yaml(id_to_smiles, args.output, args.indent)


if __name__ == '__main__':
    main()
