#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from joblib import Parallel, delayed

import polars as pl


def parse_args(nerpa_dir: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog='all_mibig_bgcs_vs_their_nrps.py',
        description=(
            'Run a single-BGC Nerpa wrapper script for every BGC that appears in pnrpdb '
            '(derived from NRP IDs like BGC0001234.5 -> BGC0001234). Each run gets its own output folder.'
        ),
    )
    parser.add_argument(
        '-o',
        '--output-dir',
        dest='output_dir',
        type=Path,
        default=nerpa_dir / 'nerpa_results/all_mibig_bgcs_vs_their_nrps',
        help='Base output directory; each BGC will be written into OUTPUT_DIR/<BGC_ID>',
    )
    parser.add_argument(
        '-t',
        '--num-threads',
        dest='num_threads',
        type=int,
        default=1,
        help='Number of parallel jobs to run (passed to joblib.Parallel; -1 uses all cores)',
    )
    return parser.parse_args()


def bgc_ids_from_pnrpdb(pnrpdb_path: Path) -> list[str]:
    # IDs look like BGC0001234.5; the corresponding BGC is the part before the dot.
    df: pl.DataFrame = pl.read_csv(pnrpdb_path, separator='\t')

    bgc_ids_df: pl.DataFrame = (
        df.filter(pl.col('ID').str.starts_with('BGC') & pl.col('ID').str.contains(r'\.'))
        .select(pl.col('ID').str.split('.').list.first().alias('bgc_id'))
        .unique()
        .sort('bgc_id')
    )

    return bgc_ids_df.get_column('bgc_id').to_list()


def main() -> None:
    nerpa_dir = Path(__file__).resolve().parent.parent
    args = parse_args(nerpa_dir)

    run_one_script: Path = (
        nerpa_dir
        / "_run_configurations"
        / "mibig_bgc_vs_its_nrps.py"
    )
    pnrpdb_path = (
        nerpa_dir
        / "data"
        / "input"
        / "pnrpdb2_expanded.tsv"
    )
    if not run_one_script.exists():
        raise FileNotFoundError(run_one_script)

    bgc_ids: list[str] = bgc_ids_from_pnrpdb(pnrpdb_path)[:5]

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f'Found {len(bgc_ids)} BGC IDs from {pnrpdb_path}. Running {run_one_script} on each one with output in {args.output_dir}...')

    def _run_one(i: int, bgc_id: str) -> tuple[str, bool]:
        out_dir: Path = args.output_dir / bgc_id
        cmd: list[str] = [
            "python",
            str(run_one_script),
            bgc_id,
            "-o",
            str(out_dir),
            "--dont-open-html",
        ]

        print(f"[{i}/{len(bgc_ids)}] {bgc_id} -> {out_dir}")
        out_dir.mkdir(parents=True, exist_ok=True)

        try:
            subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        except subprocess.CalledProcessError:
            print(f"{bgc_id} failed.", file=sys.stderr)
            return (bgc_id, False)
        else:
            print(f"{bgc_id} success.")
            return (bgc_id, True)

    Parallel(n_jobs=args.num_threads)(
        delayed(_run_one)(i, bgc_id) for i, bgc_id in enumerate(bgc_ids, start=1)
    )

    print('Done.')


if __name__ == '__main__':
    main()
