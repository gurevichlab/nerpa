import subprocess
import sys
from __future__ import annotations
from collections import Counter

import argparse
from typing import Any
from pathlib import Path
import polars as pl
from src.benchmarking.data_frames import MIBiG_BGCs_Info
from collections import defaultdict
import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Scaffold: split approved matches into train/test using MiBIG metadata, "
            "run training on train subset, then run nerpa.py on produced configs."
        )
    )
    parser.add_argument(
        "--mibig-bgcs-tsv",
        type=Path,
        required=True,
        help="Path to MiBIG BGCs TSV (metadata used to split into dissimilar train/test).",
    )
    parser.add_argument(
        "--antismash-results-dir",
        type=Path,
        required=True,
        help="Directory with antiSMASH results used by training/nerpa.",
    )
    parser.add_argument(
        "--approved-matches-yaml",
        type=Path,
        required=True,
        help="YAML file with approved matches (input set to split into train/test).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory where train/test YAMLs and other outputs will be written.",
    )
    return parser.parse_args()

import random

_SPLIT_SEED: int = 1337

def split_train_test(
    mibig_bgcs_tsv: Path,
    approved_matches_yaml: Path,
    output_dir: Path,
    split_fraction: float = 0.5,
) -> tuple[Path, Path]:
    print(f"===== Splitting approved matches into train/test with split_fraction={split_fraction}...")
    
    output_dir.mkdir(parents=True, exist_ok=True)

    with approved_matches_yaml.open("r", encoding="utf-8") as f:
        approved_matches: list[dict[str, Any]] = yaml.safe_load(f)

    mibig: MIBiG_BGCs_Info = MIBiG_BGCs_Info.from_csv(mibig_bgcs_tsv)
    # q: construct a mapping BIGSCAPE family -> list of BGCs in that family
    bgc_id_to_family = dict(
        mibig.select(
            pl.col(MIBiG_BGCs_Info.BGC_ID),
            pl.col(MIBiG_BGCs_Info.BIGSCAPE_FAMILIES),
        ).iter_rows()
    )
    bigscape_family_to_matches = defaultdict(list)
    for match in approved_matches:
        bgc_id = match["bgc_id"]
        family = bgc_id_to_family.get(bgc_id)
        if family is not None:
            bigscape_family_to_matches[family].append(match)
        else:
            print(f'BiG-SCAPE family not found for {bgc_id}')

    rng = random.Random(_SPLIT_SEED)
    fams = list(bigscape_family_to_matches.keys())
    rng.shuffle(fams)

    matches_for_training = []
    matches_for_testing = []

    total_matches = len(approved_matches)
    target_train_matches = int(total_matches * split_fraction)

    cur_cnt_for_training = 0
    fams_split_idx = 0
    for i, fam in enumerate(fams):
        if cur_cnt_for_training >= target_train_matches:
            fams_split_idx = i
            break
        cur_cnt_for_training += len(bigscape_family_to_matches[fam])

    matches_for_training = [
        match for fam in fams[:fams_split_idx]
        for match in bigscape_family_to_matches[fam]
    ]
    matches_for_testing = [
        match for fam in fams[fams_split_idx:]
        for match in bigscape_family_to_matches[fam]
    ]

    train_yaml: Path = output_dir / "approved_matches.train.yaml"
    test_yaml: Path = output_dir / "approved_matches.test.yaml"

    with train_yaml.open("w", encoding="utf-8") as f:
        yaml.safe_dump(matches_for_training, f, sort_keys=False)
    with test_yaml.open("w", encoding="utf-8") as f:
        yaml.safe_dump(matches_for_testing, f, sort_keys=False)

    print(f'Wrote {len(matches_for_training)} approved matches to {train_yaml} and {len(matches_for_testing)} to {test_yaml}.')

    return train_yaml, test_yaml


def run_training(
    nerpa_dir: Path,
    approved_matches_train_yaml: Path,
    output_dir: Path,
) -> Path:
   """
    Step 2:
    - Run another script on the training subset to generate configs / model artifacts.
    """
    train_script: Path = nerpa_dir / 'train_nerpa.py'
    training_results_dir: Path = output_dir / 'training_results'

    cmd: list[str] = [
        sys.executable,
        str(train_script),
        '--approved-matches',
        str(approved_matches_train_yaml),
        '--output-dir',
        str(training_results_dir),
    ]
    print('===== Running training:', ' '.join(cmd))
    subprocess.run(cmd, check=True)

    new_configs_dir: Path = training_results_dir / 'new_configs'
    print(f'===== Training done. Expecting new configs in: {new_configs_dir}')
    return new_configs_dir


def fetch_configs(training_output_dir: Path) -> list[Path]:
    """
    Step 3a:
    - Locate configs created by the training script.

    TODO:
    - Define what a "config" file is (YAML/JSON?) and filename pattern.
    - Decide ordering/dedup behavior.
    """
    raise NotImplementedError


def run_nerpa(
    configs: list[Path],
    antismash_results_dir: Path,
) -> None:
    """
    Step 3b:
    - Run nerpa.py for each config.

    TODO:
    - Where is nerpa.py (within repo? need an arg)?
    - Output dir / logging.
    - Parallelism or sequential.
    """
    raise NotImplementedError


def main() -> None:
    args = parse_args()
    nerpa_dir: Path = Path(__file__).resolve().parent.parent
    assert nerpa_dir.name.startswith('nerpa'), f"Expected nerpa_dir to be the nerpa repo, got {nerpa_dir}"

    # Step 1: split approved matches into training/testing subsets
    train_yaml, test_yaml = split_train_test(
        mibig_bgcs_tsv=args.mibig_bgcs_tsv,
        approved_matches_yaml=args.approved_matches_yaml,
        output_dir=args.output_dir,
    )

    # Step 2: run training on training subset
    training_output_dir = run_training(
        nerpa_dir=nerpa_dir,
        output_dir=args.output_dir,
        approved_matches_train_yaml=train_yaml,
    )
    exit(0)  # debugging

    # Step 3: fetch produced configs and run nerpa.py on them
    configs = fetch_configs(training_output_dir=training_output_dir)
    run_nerpa(configs=configs, antismash_results_dir=args.antismash_results_dir)

    # Note: test_yaml is produced but unused for now; we’ll wire it in later.


if __name__ == "__main__":
    main()
