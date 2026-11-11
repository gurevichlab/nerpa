#!/usr/bin/env python3

import argparse
import json
import sys
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nerpa-root", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--sklearn-version", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    nerpa_root: Path = args.nerpa_root.resolve()
    sys.path.insert(0, str(nerpa_root))

    import yaml

    from src.aa_specificity_prediction_model.specificity_prediction_helper import (
        SpecificityPredictionHelper,
    )
    from src.config import SpecificityPredictionConfig, load_monomer_names_helper
    from src.general_type_aliases import AA34

    with (nerpa_root / "configs" / "config.yaml").open() as config_file:
        root_config: dict[str, Any] = yaml.safe_load(config_file)

    monomer_names_helper = load_monomer_names_helper(
        nerpa_root / root_config["monomers_config"],
        nerpa_root,
    )

    specificity_config_path = (
        nerpa_root / root_config["specificity_prediction_config"]
    )
    with specificity_config_path.open() as config_file:
        specificity_config_data: dict[str, Any] = yaml.safe_load(config_file)

    specificity_config = SpecificityPredictionConfig(
        nerpa_dir=nerpa_root,
        cfg_dict=specificity_config_data,
        monomer_names_helper=monomer_names_helper,
    )
    specificity_config.paras_model = args.model.resolve()
    specificity_config.paras_model_sklearn_version = args.sklearn_version

    prediction_helper = SpecificityPredictionHelper(
        specificity_config,
        monomer_names_helper,
    )

    aa34_values: list[str] = json.load(sys.stdin)
    predictions: dict[str, dict[str, float]] = {
        aa34_value: {
            str(residue): score
            for residue, score in prediction_helper._predict_paras_a34(
                AA34(aa34_value)
            ).items()
        }
        for aa34_value in aa34_values
    }

    json.dump(predictions, sys.stdout)


if __name__ == "__main__":
    main()
