from __future__ import annotations

import argparse
from pathlib import Path


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
    return parser.parse_args()

import random

_SPLIT_SEED: int = 1337

def split_train_test(
    mibig_bgcs_tsv: Path,
    approved_matches_yaml: Path,
    output_dir: Path,
    split_fraction: float = 0.5,
) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)

    with approved_matches_yaml.open("r", encoding="utf-8") as f:
        approved_data: list[dict[str, Any]] = yaml.safe_load(f)

    approved_bgc_ids: set[str] = {item["bgc_id"] for item in approved_data}

    mibig: MIBiG_BGCs_Info = MIBiG_BGCs_Info.from_csv(mibig_bgcs_tsv)
    mibig = MIBiG_BGCs_Info(
        mibig.with_columns(
            pl.col(MIBiG_BGCs_Info.BGC_ID)
            .is_in(sorted(approved_bgc_ids))
            .alias(MIBiG_BGCs_Info.IN_APPROVED_MATCHES)
        )
    )
    approved_bgcs_info: pl.DataFrame = mibig.filter(pl.col(MIBiG_BGCs_Info.IN_APPROVED_MATCHES))

    # Count BGCs per family (families are the "dissimilarity" unit)
    fam_counts: pl.DataFrame = (
        approved_bgcs_info.group_by(MIBiG_BGCs_Info.BIGSCAPE_FAMILIES)
        .agg(pl.len().alias("n_bgcs"))
        .select([MIBiG_BGCs_Info.BIGSCAPE_FAMILIES, "n_bgcs"])
    )

    fams: list[str] = fam_counts[MIBiG_BGCs_Info.BIGSCAPE_FAMILIES].to_list()
    rng = random.Random(_SPLIT_SEED)
    rng.shuffle(fams)

    total_bgcs: int = len(approved_bgc_ids)
    target_train_bgcs: int = int(round(total_bgcs * split_fraction))
    target_train_bgcs = max(0, min(target_train_bgcs, total_bgcs))

    counts_by_fam: dict[str, int] = dict(
        zip(
            fam_counts[MIBiG_BGCs_Info.BIGSCAPE_FAMILIES].to_list(),
            fam_counts["n_bgcs"].to_list(),
            strict=True,
        )
    )

    train_fams: set[str] = set()
    running: int = 0
    for fam in fams:
        if running >= target_train_bgcs:
            break
        train_fams.add(fam)
        running += counts_by_fam.get(fam, 0)

    test_fams: set[str] = set(fams) - train_fams

    train_bgc_ids: set[str] = set(
        approved_bgcs_info.filter(pl.col(MIBiG_BGCs_Info.BIGSCAPE_FAMILIES).is_in(train_fams))
        .select(MIBiG_BGCs_Info.BGC_ID)
        .to_series()
        .to_list()
    )
    test_bgc_ids: set[str] = set(
        approved_bgcs_info.filter(pl.col(MIBiG_BGCs_Info.BIGSCAPE_FAMILIES).is_in(test_fams))
        .select(MIBiG_BGCs_Info.BGC_ID)
        .to_series()
        .to_list()
    )

    train_data: list[dict[str, Any]] = [item for item in approved_data if item["bgc_id"] in train_bgc_ids]
    test_data: list[dict[str, Any]] = [item for item in approved_data if item["bgc_id"] in test_bgc_ids]

    train_yaml: Path = output_dir / "approved_matches.train.yaml"
    test_yaml: Path = output_dir / "approved_matches.test.yaml"

    with train_yaml.open("w", encoding="utf-8") as f:
        yaml.safe_dump(train_data, f, sort_keys=False)
    with test_yaml.open("w", encoding="utf-8") as f:
        yaml.safe_dump(test_data, f, sort_keys=False)

    return train_yaml, test_yaml


def run_training(
    antismash_results_dir: Path,
    approved_matches_train_yaml: Path,
) -> Path:
    """
    Step 2:
    - Run another script on the training subset to generate configs / model artifacts.

    TODO:
    - Which training script? Provide a CLI path (maybe as an arg later).
    - Decide where artifacts are written.
    - Return path to directory containing generated configs.
    """
    raise NotImplementedError


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

    # Step 1: split approved matches into training/testing subsets
    train_yaml, test_yaml = split_train_test(
        mibig_bgcs_tsv=args.mibig_bgcs_tsv,
        approved_matches_yaml=args.approved_matches_yaml,
    )
    exit(0)  # debugging

    # Step 2: run training on training subset
    training_output_dir = run_training(
        antismash_results_dir=args.antismash_results_dir,
        approved_matches_train_yaml=train_yaml,
    )

    # Step 3: fetch produced configs and run nerpa.py on them
    configs = fetch_configs(training_output_dir=training_output_dir)
    run_nerpa(configs=configs, antismash_results_dir=args.antismash_results_dir)

    # Note: test_yaml is produced but unused for now; we’ll wire it in later.


if __name__ == "__main__":
    main()
