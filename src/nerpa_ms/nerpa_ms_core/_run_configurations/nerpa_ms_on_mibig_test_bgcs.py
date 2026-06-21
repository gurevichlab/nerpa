# ARGS="--nerpa-results /home/ilianolhin/work/tools/nerpa/nerpa_results/all_mibig_bgcs_vs_their_nrps"
import subprocess
import os
from pathlib import Path
import sys
nerpa_root = (
    Path(__file__).resolve()
    .parent.parent.parent.parent.parent
)
sys.path.append(str(nerpa_root / "src"))

from benchmarking.data_frames import (
    MIBiG_BGCs_Info,
    PNRPDB_Info,
    get_good_mibig_bgs_to_nrps,
)

import subprocess
import os
import sys
from pathlib import Path
import polars as pl
import yaml
from typing import Dict, List
import argparse
from joblib import Parallel, delayed


def parse_args(nerpa_root: Path) -> argparse.Namespace:
    nerpa_ms_root = (
        nerpa_root
        / "src"
        / "nerpa_ms"
        / "nerpa_ms_core"
    )
    parser = argparse.ArgumentParser(prog="run Nerpa on MIBiG BGC with >= min_a_domains A domains vs its NRP(s) with >= min_monomers NRPS monomers")
    parser.add_argument(
        "--nerpa-results",
        dest="nerpa_results",
        type=Path,
        help="Path to the nerpa_results directory. Each subdirectory of this directory should be a Nerpa run on a MIBiG BGC vs its NRP(s), named after the BGC ID (e.g. BGC0001234). BGCs with only one corresponding NRP are skipped",
    )
    default_output_dir = (
        nerpa_ms_root
        / "output"
        / "mibig_good_bgcs_vs_their_good_nrps"
    )
    parser.add_argument(
        "-o", "--output-dir",
        dest="output_dir",
        type=Path,
        default=default_output_dir,
        help=f"Output directory (default: {default_output_dir}"
    )

    parser.add_argument(
        "-t",
        "--num-threads",
        dest="num_threads",
        type=int,
        default=1,
        help="Number of parallel jobs to run (passed to joblib.Parallel; -1 uses all cores)",
    )

    default_max_edits = 3
    parser.add_argument(
        "--max-edits",
        dest="max_edits",
        type=int,
        default=default_max_edits,
        help=f"Maximum number of edits to consider when generating variants (default: {default_max_edits})",
    )

    default_num_variants_per_num_edits = 1000
    parser.add_argument(
        "--num-variants-per-num-edits",
        dest="num_variants_per_num_edits",
        type=int,
        default=default_num_variants_per_num_edits,
        help=f"Number of variants to generate per number of edits (default: {default_num_variants_per_num_edits})",
    )

    parser.add_argument(
        "--max-dirs-to-process",
        dest="max_dirs_to_process",
        type=int,
        default=None,
        help="Maximum number of subdirectories of --nerpa-results to process (default: unlimited). Useful for testing on a subset of the data. Used for debugging",
    )

    return parser.parse_args()

def get_nrp_ids(nerpa_results_single: Path) -> List[str]:
    nrp_representatives_path = (
        nerpa_results_single
        / "preprocessed_input"
        / "NRP_representatives.yaml"
    )
    with nrp_representatives_path.open() as f:
        nrp_representatives = yaml.safe_load(f)
    return [entry['nrp_id'][0] for entry in nrp_representatives]


def _run_one(
        nerpa_results_single: Path,
        nerpa_root: Path,
        main_output_dir: Path,
        idx: int,
        total: int,
        max_edits: int,
        num_variants_per_num_edits: int,
):
    print(f"[{idx}/{total}] Processing {nerpa_results_single.name}...")

    if len(get_nrp_ids(nerpa_results_single)) <= 1:
        print(f"Skipping {nerpa_results_single.name} because it has <= 1 NRP(s).")
        return

    bgc_id = nerpa_results_single.name
    out_dir: Path = main_output_dir / bgc_id
    args = [
        "cargo", "run", "--bin", "nerpa_ms_core", "--",
        "--nerpa-results", str(nerpa_results_single),
        "--max-edits", str(max_edits),
        "--num-variants-per-num-edits", str(num_variants_per_num_edits),
        "--out", str(out_dir),
        "--monomers-db-json", "data/monomers_db.json",
        "--nerpa-root", str(nerpa_root),
    ]
    nerpa_ms_root = (
        nerpa_root
        / 'src'
        / 'nerpa_ms'
        / 'nerpa_ms_core'
    )

    try:
        print(f"Running command:\n {' '.join(args)}")
        subprocess.run(
            args,
            cwd=nerpa_ms_root,
            check=True,
            stdout=subprocess.DEVNULL,
            #stderr=subprocess.DEVNULL,
        )
    except subprocess.CalledProcessError:
        print(f"{bgc_id} failed.", file=sys.stderr)
    else:
        print(f"{bgc_id} success.")

def main():
    os.environ["RUSTFLAGS"] = "-Awarnings"
    # os.environ["RUST_BACKTRACE"] = "1"

    args = parse_args(nerpa_root)
    single_results_dirs = []
    for results_dir in args.nerpa_results.iterdir():
        if len(get_nrp_ids(results_dir)) <= 1:
            print(f"Skipping {results_dir.name} because it has <= 1 NRP(s).")
        else:
            single_results_dirs.append(results_dir)

    if args.max_dirs_to_process is not None:
        single_results_dirs = single_results_dirs[:args.max_dirs_to_process]

    print(f"Processing {len(single_results_dirs)} directories (after filtering out those with <= 1 NRP(s)).")
    
    Parallel(n_jobs=args.num_threads)(
        delayed(_run_one)(
            nerpa_results_single=nerpa_results_one_bgc,
            nerpa_root=nerpa_root,
            main_output_dir=args.output_dir,
            idx=i,
            total=len(list(args.nerpa_results.iterdir())),
            max_edits=args.max_edits,
            num_variants_per_num_edits=args.num_variants_per_num_edits,
        )
        for i, nerpa_results_one_bgc in enumerate(single_results_dirs, start=1)
    )


if __name__ == "__main__":
    main()
