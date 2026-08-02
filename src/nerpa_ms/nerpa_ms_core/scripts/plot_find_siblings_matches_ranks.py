#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def _get_in(obj: dict[str, Any], path: list[str]) -> Any:
    cur: Any = obj
    for key in path:
        cur = cur[key]
    return cur


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            'Read a JSON list of TargetSearchResults, group found matches by '
            'linearization distance, and save histograms of matched_variant_rank per distance as SVG.'
        )
    )
    parser.add_argument(
        'json_path',
        type=Path,
        help='Path to JSON file containing a list of TargetSearchResults objects.',
    )
    parser.add_argument(
        '-o',
        '--out-dir',
        type=Path,
        default=Path('.'),
        help='Output directory for SVGs (default: current directory).',
    )
    parser.add_argument(
        '--bins',
        type=int,
        default=30,
        help='Number of histogram bins (default: 30).',
    )
    parser.add_argument(
        '--max-distance',
        type=int,
        default=None,
        help='If set, only include distances <= this value.',
    )
    args = parser.parse_args()

    # Import here so argparse --help is fast even if matplotlib isn't installed
    import matplotlib

    matplotlib.use('Agg')  # no GUI backend
    import matplotlib.pyplot as plt

    data: list[dict[str, Any]] = json.loads(args.json_path.read_text(encoding='utf-8'))

    # Collect ranks by linearization_distance, considering only entries with found_match != null
    ranks_by_dist: dict[int, list[int]] = {}
    for item in data:
        found_match = item.get('found_match')
        if found_match is None:
            continue

        dist_any = _get_in(item, ['features', 'linearization_distance'])
        dist: int = int(dist_any)
        if args.max_distance is not None and dist > args.max_distance:
            continue

        rank_any = found_match.get('matched_variant_rank')
        if rank_any is None:
            continue
        rank: int = int(rank_any)

        ranks_by_dist.setdefault(dist, []).append(rank)

    if not ranks_by_dist:
        raise SystemExit('No found matches present (or no usable matched_variant_rank fields).')

    args.out_dir.mkdir(parents=True, exist_ok=True)

    for dist in sorted(ranks_by_dist.keys()):
        ranks = ranks_by_dist[dist]

        plt.figure(figsize=(8, 4.5))
        plt.hist(ranks, bins=args.bins)
        plt.title(f'matched_variant_rank histogram (linearization_distance={dist}, n={len(ranks)})')
        plt.xlabel('matched_variant_rank')
        plt.ylabel('count')
        plt.tight_layout()

        out_path = args.out_dir / f'matches_ranks_hist_dist{dist}.svg'
        plt.savefig(out_path, format='svg')
        plt.close()

    print(f'Saved {len(ranks_by_dist)} SVG(s) to {args.out_dir}')


if __name__ == '__main__':
    main()
