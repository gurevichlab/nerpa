from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

import yaml

PARAMS: dict[str, Any] = {
    # Paths relative to nerpa_dir (script_location/..)
    "approved_matches_rel": "data/for_training_and_testing/approved_matches.yaml",
    "mibig_bgcs_info_rel": "data/for_training_and_testing/mibig_bgcs_info.tsv",
    "local_paths_rel": "local_paths.yaml",
    "output_dir_rel": "benchmarking/alignment_reconstruction",
    # Key inside local_paths.yaml
    "antismash_results_key": "antismash_results",
    # Other CLI params
    "split_fraction": 0.5,
    # Main script filename (assumed next to this helper script)
    "benchmark_script": "benchmark_alignment_reconstruction_main.py",
}


def main() -> None:
    script_dir: Path = Path(__file__).resolve().parent
    nerpa_dir: Path = script_dir / ".."

    approved_matches: Path = nerpa_dir / PARAMS["approved_matches_rel"]
    mibig_bgcs_info: Path = nerpa_dir / PARAMS["mibig_bgcs_info_rel"]
    output_dir: Path = nerpa_dir / PARAMS["output_dir_rel"]

    local_paths_path: Path = nerpa_dir / PARAMS["local_paths_rel"]
    with local_paths_path.open("r", encoding="utf-8") as f:
        local_paths: dict[str, Any] = yaml.safe_load(f)

    antismash_results = Path(local_paths[PARAMS["antismash_results_key"]])

    benchmark_script: Path = script_dir / PARAMS["benchmark_script"]

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

    subprocess.run(["python", str(benchmark_script), *args], check=True)
