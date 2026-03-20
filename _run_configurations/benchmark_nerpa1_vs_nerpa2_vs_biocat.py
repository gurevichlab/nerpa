from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any
import os
import subprocess
import sys

import yaml
import argparse

PARAMS: dict[str, Any] = {
    # Paths relative to nerpa_dir (script_location/..)
    "nerpa1_results": "data/for_training_and_testing/nerpa1_report_mibig4_vs_mibig_norine.csv",
    "nerpa2_results": "nerpa_results/mibig4_vs_mibig_norine",
    "biocat_results": "data/for_training_and_testing/biocat_results.txt",

    "local_paths_rel": "local_paths.yaml",
    "output_dir_rel": "benchmarking",

    # Key inside local_paths.yaml
    "antismash_results_key": "as_results_mibig4_nrps",

    # Main script filename
    "benchmark_script": "benchmarking_plots.py",
}


def run_nerpa2(nerpa_exec: Path,
               smiles_tsv: Path,
               antismash_results: Path,
               output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    command = [
        sys.executable, str(nerpa_exec),
        "--antismash", str(antismash_results),
        "--smiles-tsv", str(smiles_tsv),
        "--col-id", "ID",
        "--output-dir", str(output_dir),
        "--force-output-dir",
        "--fast-matching",
        '--min-num-matches-per-bgc', '10',
        '--max-num-matches-per-bgc', '10',
        '--min-num-matches-per-nrp', '0',
        '--max-num-matches-per-nrp', '0',
        '--max-num-matches', '0',
        '--disable-bgc-deduplication',
        '--process-hybrids',
        '--threads', '8',
        '--skip-molecule-drawing',
#        "--let-it-crash",
        "--dump-all-preprocessed",
    ]
    # print(f"Running Nerpa 2 with command: {' '.join(command)}")
    subprocess.run(command, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark Nerpa 1 vs Nerpa 2 vs Biocat on MiBIG4 vs MiBIG-Norine")
    parser.add_argument('--rerun-nerpa2', action='store_true', help="Whether to rerun Nerpa 2 or use existing results")
    return parser.parse_args()

def main() -> None:
    script_dir: Path = Path(__file__).resolve().parent
    nerpa_dir: Path = script_dir.parent
    print(f"Identified nerpa_dir={nerpa_dir}")
    args = parse_args()

    with open(nerpa_dir / PARAMS["local_paths_rel"]) as f:
        local_paths = yaml.safe_load(f)

    antismash_results = Path(local_paths[PARAMS["antismash_results_key"]])

    if args.rerun_nerpa2:
        print("Rerunning Nerpa 2 on MiBIG4 vs MiBIG-Norine...")
        run_nerpa2(nerpa_exec=nerpa_dir / "nerpa.py",
                   smiles_tsv=nerpa_dir / "data/input/pnrpdb2_mibig_norine_at_least3_mons.tsv",
                   antismash_results=antismash_results,
                   output_dir=nerpa_dir / PARAMS["nerpa2_results"])
    else:
        print(f'Using existing Nerpa 2 results in {nerpa_dir / PARAMS["nerpa2_results"]}')


    benchmark_script: Path = nerpa_dir / PARAMS["benchmark_script"]

    command_args: list[str] = [
        '--nerpa1-report', str(nerpa_dir / PARAMS["nerpa1_results"]),
        '--nerpa2-report', str(nerpa_dir / PARAMS["nerpa2_results"] / "report.tsv"),
        '--biocat-report', str(nerpa_dir / PARAMS["biocat_results"]),
        "--output-dir", str(nerpa_dir / PARAMS["output_dir_rel"]),
        #"--split-fraction",
        #str(PARAMS["split_fraction"]),
    ]

    env: dict[str, str] = dict(os.environ)
    env["PYTHONPATH"] = str(nerpa_dir)

    command: list[str] = [sys.executable, str(benchmark_script), *command_args]
    print(f"Running command: PYTHONPATH={env['PYTHONPATH']} {sys.executable} {benchmark_script} {' '.join(command_args)}")

    subprocess.run(command, check=True, env=env)

if __name__ == "__main__":
    main()
