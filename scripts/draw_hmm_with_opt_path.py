#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Optional
from src.hmm.hmm import HMM


def _parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description='Draw an HMM and optionally highlight an optimal path.')
    p.add_argument('--input', type=Path, required=True, help='Path to JSON produced by Rust code')
    p.add_argument('--output', type=Path, required=True, help='Output image filename (passed to HMM.draw)')
    return p.parse_args(argv)


def main(argv: list[str]) -> int:
    args = _parse_args(argv)

    input_path: Path = args.input
    output_path: Path = args.output

    with input_path.open('r', encoding='utf-8') as f:
        payload: dict[str, Any] = json.load(f)

    hmm_dict: dict[str, Any] = payload['hmm']
    opt_path_raw: list[Any] = payload.get('opt_path', [])
    emission_labels: Optional[list[str]] = payload.get('emission_labels')

    opt_path: list[int] = [int(x) for x in opt_path_raw]

    hmm = HMM.from_json_dict(hmm_dict)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    hmm.draw(
        filename=output_path,
        highlight_path=opt_path if opt_path else None,
        emission_names=emission_labels if emission_labels else None,
    )

    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1:]))
