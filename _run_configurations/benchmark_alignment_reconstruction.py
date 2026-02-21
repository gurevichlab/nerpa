from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any
import os
import subprocess
import sys

import yaml

PARAMS: dict[str, Any] = {
    # Paths relative to nerpa_dir (script_location/..)
    "approved_matches_rel": "data/for_training_and_testing/approved_matches.yaml",
    "mibig_bgcs_info_rel": "data/for_training_and_testing/mibig_bgcs_info.tsv",
    "local_paths_rel": "local_paths.yaml",
    "output_dir_rel": "benchmarking/alignment_reconstruction",
    # Key inside local_paths.yaml
    "antismash_results_key": "as_results_mibig4_nrps",
    # Other CLI params
    "num_groups": 5,
    # Main script filename (assumed next to this helper script)
    "benchmark_script": "src/benchmarking/benchmark_alignment_reconstruction_main.py",
}


def main() -> None:
    script_dir: Path = Path(__file__).resolve().parent
    nerpa_dir: Path = script_dir.parent
    print(f"Identified nerpa_dir={nerpa_dir}")

    approved_matches: Path = nerpa_dir / PARAMS["approved_matches_rel"]
    mibig_bgcs_info: Path = nerpa_dir / PARAMS["mibig_bgcs_info_rel"]
    output_dir: Path = nerpa_dir / PARAMS["output_dir_rel"]

    local_paths_path: Path = nerpa_dir / PARAMS["local_paths_rel"]
    with local_paths_path.open("r", encoding="utf-8") as f:
        local_paths: dict[str, Any] = yaml.safe_load(f)

    antismash_results = Path(local_paths[PARAMS["antismash_results_key"]])

    benchmark_script: Path = nerpa_dir / PARAMS["benchmark_script"]

    args: list[str] = [
        "--mibig-bgcs-tsv",
        str(mibig_bgcs_info),
        "--antismash-results-dir",
        str(antismash_results),
        "--approved-matches-yaml",
        str(approved_matches),
        "--output-dir",
        str(output_dir),
        #"--split-fraction",
        #str(PARAMS["split_fraction"]),
    ]

    env: dict[str, str] = dict(os.environ)
    env["PYTHONPATH"] = str(nerpa_dir)

    command: list[str] = [sys.executable, str(benchmark_script), *args]
    print(f"Running command: PYTHONPATH={env['PYTHONPATH']} {sys.executable} {benchmark_script} {' '.join(args)}")

    subprocess.run(command, check=True, env=env)

if __name__ == "__main__":
    main()
