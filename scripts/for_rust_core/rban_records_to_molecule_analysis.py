#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, TypedDict

import yaml

from src.build_output.chem_helper import MolRecord
from src.rban_parsing.rban_parser import Parsed_rBAN_Record


class CanonicalSmilesData(TypedDict):
    smiles: str
    atom_ids_in_smiles_order: list[int]


def build_id_to_canonical_smiles_data(
    records_raw: list[Any],
) -> dict[str, CanonicalSmilesData]:
    id_to_data: dict[str, CanonicalSmilesData] = {}

    for record_raw in records_raw:
        record = Parsed_rBAN_Record.from_dict(record_raw)
        mol_record = MolRecord.from_rban_record(record)

        id_to_data[str(record.compound_id)] = {
            "smiles": mol_record.to_canonical_smiles(),
            "atom_ids_in_smiles_order": mol_record.canonical_atom_order(),
        }

    return id_to_data


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Analyze Parsed_rBAN_Record objects with RDKit and return canonical "
            "SMILES and canonical atom order."
        )
    )
    parser.add_argument(
        "input",
        type=Path,
        help="YAML file containing a list of Parsed_rBAN_Record objects.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Write output YAML here instead of stdout.",
    )
    return parser.parse_args()


def main() -> None:
    import sys

    args = parse_args()
    records_raw: Any = yaml.safe_load(
        args.input.read_text(encoding="utf-8")
    )
    if not isinstance(records_raw, list):
        raise SystemExit("Input YAML must be a list.")

    output = yaml.safe_dump(
        build_id_to_canonical_smiles_data(records_raw),
        allow_unicode=True,
        sort_keys=False,
    )

    if args.output is None:
        sys.stdout.write(output)
    else:
        args.output.write_text(output, encoding="utf-8")


if __name__ == "__main__":
    main()