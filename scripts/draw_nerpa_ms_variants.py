#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import List, Optional, Tuple, Any, NamedTuple
from src.rban_parsing.rban_parser import (
    Parsed_rBAN_Record,
    MonomerIdx,
)
from src.general_type_aliases import LogProb
from collections import defaultdict
from src.config import load_monomer_names_helper
from src.build_output.draw_graph import (
    draw_molecule_diff,
    draw_monomer_graph_diff,
)

def _parse_opt_int(x: Any) -> Optional[int]:
    if x is None:
        return None
    if isinstance(x, bool):
        raise ValueError(f"Expected int or None, got bool: {x!r}")
    if isinstance(x, int):
        return x
    if isinstance(x, float):
        if x.is_integer():
            return int(x)
        raise ValueError(f"Expected integer-valued float or None, got: {x!r}")
    if isinstance(x, str):
        s = x.strip()
        if s == "" or s.lower() in {"none", "null"}:
            return None
        return int(s)
    raise ValueError(f"Expected int/str/None, got: {type(x).__name__}: {x!r}")


class Altered_rBAN_Record(NamedTuple):
    score: LogProb
    new_record: Parsed_rBAN_Record
    old_to_new_mon_map: List[Tuple[Optional[MonomerIdx], Optional[MonomerIdx]]]

    @classmethod
    def from_dict(cls, d: dict) -> Altered_rBAN_Record:
        score = float(d["score"])
        new_record = Parsed_rBAN_Record.from_dict(d["new_record"])

        raw_map = d.get("old_to_new_mon_map", [])
        old_to_new_mon_map: List[Tuple[Optional[int], Optional[int]]] = [
            (_parse_opt_int(pair[0]), _parse_opt_int(pair[1]))
            for pair in raw_map
        ]

        return cls(
            score=score,
            new_record=new_record,
            old_to_new_mon_map=old_to_new_mon_map,
        )

class ItemForDrawing(NamedTuple):
    original: Parsed_rBAN_Record
    new_variant: Altered_rBAN_Record

    @classmethod
    def from_dict(cls, d: dict) -> ItemForDrawing:
        original = Parsed_rBAN_Record.from_dict(d["original"])
        new_variant = Altered_rBAN_Record.from_dict(d["new_variant"])
        return cls(original=original, new_variant=new_variant)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_json",
                        required=True,
                        type=Path,
                        help="Path to items_for_drawing.json")
    parser.add_argument("--output_dir",
                        required=True,
                        type=Path,
                        help="Directory to write drawings")
    return parser.parse_args()

def main() -> None:
    nerpa_dir = Path(__file__).parent.parent.parent
    assert (nerpa_dir / "nerpa.py").exists(), f"Invalid nerpa_dir: {nerpa_dir}"
    monomers_cfg_file = nerpa_dir / "configs/monomers_config.yaml"
    monomer_names_helper = load_monomer_names_helper(monomers_cfg_file, nerpa_dir)

    args = parse_arguments()


    data = json.loads(args.input_json.read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError("Expected input JSON to be a list")
    items_for_drawing = [ItemForDrawing.from_dict(item) for item in data]

    items_by_compound_id = defaultdict(list)
    for item in items_for_drawing:
        items_by_compound_id[item.original.compound_id].append(item)

    for compound_id, items in items_by_compound_id.items():
        out_dir = args.output_dir / compound_id
        out_dir.mkdir(parents=True, exist_ok=True)

        for i, item in enumerate(items):
            draw_monomer_graph_diff(
                original=item.original,
                modified=item.new_variant.new_record,
                old_to_new_mon_map=item.new_variant.old_to_new_mon_map,
                monomer_names_helper=monomer_names_helper,
                output=out_dir / f"monomer_graph_diff_{compound_id}_{i:03d}.svg",)
            draw_molecule_diff(
                original=item.original,
                modified=item.new_variant.new_record,
                old_to_new_mon_map=item.new_variant.old_to_new_mon_map,
                monomer_names_helper=monomer_names_helper,
                output=out_dir / f"monomer_graph_diff_{compound_id}_{i:03d}.svg",)


if __name__ == "__main__":
    main()
