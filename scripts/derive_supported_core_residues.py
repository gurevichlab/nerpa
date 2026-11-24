#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import Any

import polars as pl
import yaml

from src.monomer_names_helper.monomer_names_helper import MonomerNamesHelper
from src.monomer_names_helper.monomer_names_table import MonomerNameFormat
from src.paras.paras_wrapper import ParasWrapper


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Derive Nerpa's supported core-residue set."
    )
    parser.add_argument(
        "--nerpa-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Nerpa repository root.",
    )
    return parser.parse_args()


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open() as file:
        return yaml.safe_load(file)


def main() -> None:
    args = parse_args()
    nerpa_dir: Path = args.nerpa_dir.resolve()

    specificity_config = load_yaml(
        nerpa_dir / "configs/specificity_prediction_config.yaml"
    )
    monomers_config = load_yaml(nerpa_dir / "configs/monomers_config.yaml")

    model = ParasWrapper(
        model_dump=nerpa_dir / specificity_config["paras_model"],
        sklearn_version=specificity_config["paras_model_sklearn_version"],
    )
    paras_predictions_raw = model.predict("A" * 34)
    paras_residues = sorted(map(str, paras_predictions_raw))

    print(f"PARAS substrates:\n {paras_residues}")
    print(f"Number of raw PARAS substrates: {len(paras_residues)}")

    monomer_names_helper = MonomerNamesHelper(
        names_table_path=nerpa_dir / monomers_config["monomer_names_table"],
        supported_residues=None,
    )

    paras_to_nerpa_residue: dict[str, str] = {
        residue: str(
            monomer_names_helper.parsed_name(
                residue,
                name_format=MonomerNameFormat.PARAS,
            ).residue
        )
        for residue in paras_residues
    }
    paras_nerpa_residues = sorted(set(paras_to_nerpa_residue.values()))

    print(
        "\nPARAS substrates converted to Nerpa core residues:"
        f"\n {paras_nerpa_residues}"
    )
    print(
        "Number of distinct Nerpa core residues among PARAS substrates: "
        f"{len(paras_nerpa_residues)}"
    )

    substrates_by_core_residue: dict[str, list[str]] = {}
    for substrate, core_residue in paras_to_nerpa_residue.items():
        substrates_by_core_residue.setdefault(core_residue, []).append(substrate)

    collapsed_substrates = {
        core_residue: substrates
        for core_residue, substrates in substrates_by_core_residue.items()
        if len(substrates) > 1
    }

    print("\nPARAS substrates collapsed to the same Nerpa core residue:")
    for core_residue, substrates in sorted(collapsed_substrates.items()):
        print(f"{core_residue}: {', '.join(sorted(substrates))}")

    norine_monomers_path = (
        nerpa_dir / monomers_config["norine_monomers_info"]
    )
    norine_residue_frequencies = (
        pl.read_csv(
            norine_monomers_path,
            separator="\t",
            columns=["Nerpa_residue", "count"],
        )
        .drop_nulls("Nerpa_residue")
        .group_by("Nerpa_residue")
        .agg(pl.col("count").sum())
        .with_columns(
            (pl.col("count") / pl.col("count").sum()).alias("frequency")
        )
        .sort("Nerpa_residue")
    )

    norine_residues: list[str] = norine_residue_frequencies[
        "Nerpa_residue"
    ].to_list()
    print(f"\nNorine database contains {len(set(norine_residues))} distinct core residues.")

    unsupported_paras_residues = sorted(
        set(paras_nerpa_residues) - set(norine_residues)
    )
    print("\nPARAS residues absent from the Norine database:")
    print(unsupported_paras_residues)
    print(
        "Number of PARAS residues absent from the Norine database: "
        f"{len(unsupported_paras_residues)}"
    )

    supported_core_residues = sorted(
        set(paras_nerpa_residues) & set(norine_residues)
    )
    print(f"\nSupported core residues: {supported_core_residues}")
    print(
        f"Number of supported core residues: {len(supported_core_residues)}"
    )


if __name__ == "__main__":
    main()
