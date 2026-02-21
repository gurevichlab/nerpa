from __future__ import annotations

import math
import subprocess
import sys
from collections import Counter

import argparse
from itertools import chain
from typing import Any, List
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
    parser.add_argument("--num-groups",
                        type=int,
                        default=5,
                        help="Number of groups to split into (default: 5)")
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
    num_groups: int = 5,
) -> List[List[dict]]:
    print(f"===== Splitting approved matches into {num_groups} equal sized dissimilar...")
    
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

    total_matches = len(approved_matches)
    group_size = math.ceil(total_matches / num_groups)
    groups = [[] for _ in range(num_groups)]
    cur_group_idx = 0

    for fam in fams:
        groups[cur_group_idx].extend(bigscape_family_to_matches[fam])
        if len(groups[cur_group_idx]) >= group_size:
            cur_group_idx += 1

    print('Group sizes:', [len(group) for group in groups])  # print group sizes for debugging

    # group_yamls = [
    #     output_dir / f"approved_matches.group{idx+1}.yaml"
    #     for idx in range(num_groups)
    # ]
    # for group, group_yaml in zip(groups, group_yamls):
    #     with group_yaml.open("w", encoding="utf-8") as f:
    #         yaml.safe_dump(group, f, sort_keys=False)

    return groups


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


def run_nerpa(
        nerpa_dir: Path,
        configs: Path,
        bgc_ids: set[str],
        smiles_tsv: Path,
        antismash_results_dir: Path,
        output_dir: Path,
) -> Path:
    # 1. Prepare antismash-paths-file
    antismash_paths_file: Path = output_dir / 'antismash_paths.txt'
    with antismash_paths_file.open('w', encoding='utf-8') as f:
        for bgc_id in bgc_ids:
            antismash_result_path = antismash_results_dir / bgc_id
            if not antismash_result_path.exists():
                print(f'Warning: antiSMASH result for {bgc_id} not found at {antismash_result_path}')
                continue
            f.write(str(antismash_result_path) + '\n')

    nerpa_script = nerpa_dir / 'nerpa.py'
    cmd: list[str] = [
        sys.executable, str(nerpa_script),
        '--antismash-paths-file', str(antismash_paths_file),
        '--smiles-tsv', str(smiles_tsv),
        '--col-id', 'ID',
        '--output-dir', str(output_dir / 'nerpa_results'),
        '--configs-dir', str(configs),
        '--min-num-matches-per-bgc', '50',
        '--max-num-matches-per-bgc', '50',
        '--max-num-matches', '0',
        '--force-output-dir',
        '--fast-matching',
        '--let-it-crash',
        '--dump-all-preprocessed',
        '--skip-molecule-drawing',
        '--disable-deduplication',
        '--threads', '8',
    ]

    print('===== Running nerpa.py on test set with configs from training:', ' '.join(cmd))
    subprocess.run(cmd, check=True)

    return output_dir / 'nerpa_results'


def benchmark_nerpa1_vs_2(
    nerpa_dir: Path,
    approved_matches_yaml: Path,
    nerpa2_matches: Path,
    plots_dir: Path,
) -> Path:
    nerpa1_results = (
            nerpa_dir
            / 'data'
            / 'for_training_and_testing'
            / 'nerpa1_results'
    )
    benchmark_script = (
            nerpa_dir
            / 'src'
            / 'benchmarking'
            / 'benchmark_alignment_reconstruction.py'
    )

    cmd: list[str] = [
        sys.executable, str(benchmark_script),
        '--approved-matches-yaml', str(approved_matches_yaml),
        '--nerpa1-results', str(nerpa1_results),
        '--nerpa2-matches', str(nerpa2_matches),
        '--plots-dir', str(plots_dir),
    ]
    print('===== Benchmarking Nerpa 1 vs Nerpa 2 on test set:', ' '.join(cmd))
    subprocess.run(cmd, check=True)

    return plots_dir



def main() -> None:
    args = parse_args()
    nerpa_dir: Path = Path(__file__).resolve().parent.parent.parent
    assert nerpa_dir.name.startswith('nerpa'), f"Expected nerpa_dir to be the nerpa repo, got {nerpa_dir}"

    # Step 1: split approved matches into training/testing subsets
    groups = split_train_test(
        mibig_bgcs_tsv=args.mibig_bgcs_tsv,
        approved_matches_yaml=args.approved_matches_yaml,
        output_dir=args.output_dir,
    )

    # Step 2: run training on training subset
    for i, test_group in enumerate(groups):
        print(f"===== Processing group {i+1}/{len(groups)}: {len(test_group)} matches in test set, {sum(len(g) for j,g in enumerate(groups) if j!=i)} matches in train set")
        train_group = list(chain(
            *[groups[j] for j in range(len(groups)) if j != i]
        ))

        test_yaml = args.output_dir / f"approved_matches_test_{i+1}.yaml"
        train_yaml = args.output_dir / f"approved_matches_train_{i+1}.yaml"
        with test_yaml.open("w", encoding="utf-8") as f:
            yaml.safe_dump(test_group, f, sort_keys=False)
        with train_yaml.open("w", encoding="utf-8") as f:
            yaml.safe_dump(list(train_group), f, sort_keys=False)

        print("Train/test YAMLs written. Running training on train set...")
        new_configs = run_training(
            nerpa_dir=nerpa_dir,
            output_dir=args.output_dir / f"{i+1}",
            approved_matches_train_yaml=train_yaml,
        )
        print(f"Training completed. New configs should be in {new_configs}. Now running nerpa.py on test set with these configs...")

        # Step 3: run nerpa.py on test subset using configs/artifacts from training
        bgc_ids = set(match["bgc_id"]
                      for match in yaml.safe_load(test_yaml.read_text(encoding="utf-8")))
        run_nerpa(nerpa_dir=nerpa_dir,
                  configs=new_configs,
                  bgc_ids=bgc_ids,
                  smiles_tsv=nerpa_dir / 'data/input/pnrpdb2_mibig_norine.tsv',
                  antismash_results_dir=args.antismash_results_dir,
                  output_dir=args.output_dir / f"{i+1}",)
        print("Done.")

    # q: collect Nerpa 2 matches from all test groups into a single file for benchmarking
    print('Collecting Nerpa 2 matches from all test groups into a single file for benchmarking...')
    nerpa2_matches = []
    for i in range(len(groups)):
        nerpa2_matches_file = (
                args.output_dir
                / f"{i+1}"
                / 'nerpa_results'
                / 'matches_details'
                / 'matches.yaml'
        )
        with nerpa2_matches_file.open("r", encoding="utf-8") as f:
            new_matches = yaml.safe_load(f)
            nerpa2_matches.extend(new_matches)

    with (args.output_dir / 'nerpa2_joined_matches.yaml').open("w", encoding="utf-8") as f:
        yaml.safe_dump(nerpa2_matches, f, sort_keys=False)

    print('Nerpa 2 matches collected and written to:', args.output_dir / 'nerpa2_joined_matches.yaml')

    print("All groups processed. Now benchmarking Nerpa 1 vs Nerpa 2 on test set...")
    # Step 4: benchmark Nerpa 1 vs Nerpa 2 performance on test set
    benchmark_nerpa1_vs_2(
        nerpa_dir=nerpa_dir,
        approved_matches_yaml=args.approved_matches_yaml,
        nerpa2_matches=args.output_dir / 'nerpa2_joined_matches.yaml',
        plots_dir=args.output_dir / 'plots',
    )
    print("Done.")


if __name__ == "__main__":
    main()
