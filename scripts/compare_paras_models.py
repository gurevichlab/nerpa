#!/usr/bin/env python3

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, TypedDict

import matplotlib.pyplot as plt
import polars as pl
import yaml
from sklearn.metrics import average_precision_score, precision_recall_curve
from src.config import SpecificityPredictionConfig, load_monomer_names_helper
from src.hmm.hmm_constructor.hmm_constructor_state_edge_context_relations import (
    MATCHING_STATE_TYPES,
)
from src.monomer_names_helper import (
    UNKNOWN_RESIDUE,
    MonomerNamesHelper,
    NerpaResidue,
    NRP_Monomer,
)
from src.aa_specificity_prediction_model.specificity_prediction_helper import (
    AA34,
    SpecificityPredictionHelper,
)


class ModelRunningSpec(TypedDict):
    model_path: Path
    python_path: Path
    sklearn_version: str


class ModelEvaluation(TypedDict):
    model_name: str
    labels: list[bool]
    scores: list[float]
    average_precision: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare two PARAS models using precision-recall curves. "
            "Each model can use a different Python environment."
        )
    )
    parser.add_argument(
        "--table",
        type=Path,
        required=True,
        help="TSV with EMISSION, AA34, and STATE_TYPE columns",
    )
    parser.add_argument(
        "--model-1",
        type=Path,
        required=True,
        help="Path to the first PARAS model",
    )
    parser.add_argument(
        "--python-1",
        type=Path,
        required=True,
        help="Python interpreter compatible with the first model",
    )
    parser.add_argument(
        "--sklearn-version-1",
        required=True,
        help="sklearn version compatible with the first model",
    )
    parser.add_argument(
        "--model-2",
        type=Path,
        required=True,
        help="Path to the second PARAS model",
    )
    parser.add_argument(
        "--python-2",
        type=Path,
        required=True,
        help="Python interpreter compatible with the second model",
    )
    parser.add_argument(
        "--sklearn-version-2",
        required=True,
        help="sklearn version compatible with the second model",
    )
    parser.add_argument(
        "--nerpa-root",
        type=Path,
        required=True,
        help="Path to the NERPA repository root",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Output directory",
    )
    return parser.parse_args()


def load_true_specificities(table_path: Path) -> dict[AA34, NerpaResidue]:
    matching_state_names = [
        state_type.name for state_type in MATCHING_STATE_TYPES
    ]

    # Only matching HMM states contribute amino-acid specificity observations.
    table = (
        pl.read_csv(table_path, separator="\t")
        .filter(
            pl.col("EMISSION").is_not_null()
            & pl.col("AA34").is_not_null()
            & (pl.col("EMISSION") != "")
            & (pl.col("AA34") != "")
            & pl.col("STATE_TYPE").is_in(matching_state_names)
        )
        .select("EMISSION", "AA34")
    )

    return {
        AA34(aa34): NRP_Monomer.from_str_short(str(emission)).residue
        for emission, aa34 in table.iter_rows()
    }


def load_prediction_helper(
        nerpa_root: Path,
        model_path: Path,
        sklearn_version: str,
) -> SpecificityPredictionHelper:
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

    # substitute the model path and sklearn version for this prediction helper
    specificity_config.paras_model = model_path.resolve()
    specificity_config.paras_model_sklearn_version = sklearn_version

    return SpecificityPredictionHelper(
        specificity_config,
        monomer_names_helper,
    )


def get_paras_predictions(
    aa34s: list[AA34],
    running_spec: ModelRunningSpec,
    nerpa_root: Path,
) -> dict[AA34, dict[NerpaResidue, float]]:
    worker_path = nerpa_root / "scripts" / "predict_paras_model.py"

    result = subprocess.run(
        [
            str(running_spec["python_path"]),
            str(worker_path),
            "--nerpa-root",
            str(nerpa_root),
            "--model",
            str(running_spec["model_path"]),
            "--sklearn-version",
            running_spec["sklearn_version"],
        ],
        input=json.dumps([str(aa34) for aa34 in aa34s]),
        stdout=subprocess.PIPE,
        text=True,
        check=True,
    )

    raw_predictions: dict[str, dict[str, float]] = json.loads(result.stdout)

    return {
        AA34(aa34): {
            NerpaResidue(residue): score
            for residue, score in residue_scores.items()
        }
        for aa34, residue_scores in raw_predictions.items()
    }
    
def collect_points(
    true_specs: dict[AA34, NerpaResidue],
    predictions: dict[AA34, dict[NerpaResidue, float]],
) -> tuple[list[bool], list[float]]:
    labels: list[bool] = []
    scores: list[float] = []

    for aa34, actual_residue in true_specs.items():
        # Every residue score is an independent binary observation.
        for predicted_residue, score in predictions[aa34].items():
            labels.append(predicted_residue == actual_residue)
            scores.append(score)

    return labels, scores


def plot_precision_recall(
    evaluations: list[ModelEvaluation],
    output_path: Path,
) -> None:
    figure, axis = plt.subplots(figsize=(7, 6))

    for evaluation in evaluations:
        precision, recall, _ = precision_recall_curve(
            evaluation["labels"],
            evaluation["scores"],
        )
        axis.plot(
            recall,
            precision,
            linewidth=2,
            label=(
                f"{evaluation['model_name']} "
                f"(AP={evaluation['average_precision']:.3f})"
            ),
        )

    positive_fraction = (
        sum(evaluations[0]["labels"]) / len(evaluations[0]["labels"])
    )
    axis.axhline(
        positive_fraction,
        color="gray",
        linestyle="--",
        linewidth=1,
        label=f"Random baseline ({positive_fraction:.3f})",
    )

    axis.set(
        xlabel="Recall",
        ylabel="Precision",
        title="PARAS model precision-recall comparison",
        xlim=(0.0, 1.0),
        ylim=(0.0, 1.02),
    )
    axis.grid(alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(output_path, dpi=300)
    plt.close(figure)


def plot_score_correctness(
    evaluations: list[ModelEvaluation],
    output_path: Path,
    bin_count: int = 10,
    min_points_per_bin: int = 10,
) -> None:
    bin_width = 1.0 / bin_count
    bin_centers = [
        bin_index * bin_width + bin_width / 2
        for bin_index in range(bin_count)
    ]
    bar_width = bin_width * 0.8 / len(evaluations)

    figure, axis = plt.subplots(figsize=(11, 6))

    for model_index, evaluation in enumerate(evaluations):
        bin_totals: list[int] = [0] * bin_count
        bin_correct: list[int] = [0] * bin_count

        # Group binary correctness observations by their predicted score.
        for label, score in zip(
            evaluation["labels"],
            evaluation["scores"],
            strict=True,
        ):
            bin_index = min(max(int(score / bin_width), 0), bin_count - 1)
            bin_totals[bin_index] += 1
            bin_correct[bin_index] += int(label)

        offset = (
            model_index - (len(evaluations) - 1) / 2
        ) * bar_width
        populated_bins = [
            (center, correct / total)
            for center, correct, total in zip(
                bin_centers,
                bin_correct,
                bin_totals,
                strict=True,
            )
            if total >= min_points_per_bin
        ]

        axis.bar(
            [center + offset for center, _ in populated_bins],
            [correctness for _, correctness in populated_bins],
            width=bar_width,
            label=evaluation["model_name"],
        )

    axis.set(
        xlabel="Prediction score bin",
        ylabel="Fraction correct",
        title="Prediction correctness by score",
        xlim=(0.0, 1.0),
        ylim=(0.0, 1.02),
    )
    axis.set_xticks(
        bin_centers,
        [
            f"{bin_index * bin_width:.2f}–{(bin_index + 1) * bin_width:.2f}"
            for bin_index in range(bin_count)
        ],
        rotation=90,
    )

    axis.plot(
        [0.0, 1.0],
        [0.0, 1.0],
        color="gray",
        linestyle="--",
        linewidth=1,
        label="Perfect calibration",
    )

    axis.grid(axis="y", alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(output_path, dpi=300)
    plt.close(figure)


def main() -> None:
    args = parse_args()
    nerpa_root = args.nerpa_root.resolve()
    sys.path.insert(0, str(nerpa_root))

    true_specs: dict[AA34, NerpaResidue] = load_true_specificities(args.table)


    model_running_specs: list[ModelRunningSpec] = [
        {
            "model_path": args.model_1,
            "python_path": args.python_1,
            "sklearn_version": args.sklearn_version_1,
        },
        {
            "model_path": args.model_2,
            "python_path": args.python_2,
            "sklearn_version": args.sklearn_version_2,
        },
    ]

    evaluations: list[ModelEvaluation] = []

    for running_spec in model_running_specs:
        predictions = get_paras_predictions(
            list(true_specs.keys()),
            running_spec,
            nerpa_root,
        )
        labels, scores = collect_points(true_specs, predictions)
        average_precision = average_precision_score(labels, scores)

        evaluations.append(
            {
                "model_name": running_spec["model_path"].name,
                "labels": labels,
                "scores": scores,
                "average_precision": average_precision,
            }
        )

        top_1_accuracy = sum(
            max(
                predictions[aa34],
                key=lambda residue: predictions[aa34][residue],
            )
            == actual_residue
            for aa34, actual_residue in true_specs.items()
        ) / len(true_specs)
        mean_true_score = sum(
            predictions[aa34][actual_residue]
            for aa34, actual_residue in true_specs.items()
        ) / len(true_specs)

        print(f"\n{running_spec['model_path']}")
        print(f"  Unique signatures: {len(true_specs)}")
        print(f"  Average precision: {average_precision:.4f}")
        print(f"  Top-1 accuracy:    {top_1_accuracy:.4f}")
        print(f"  Mean true score:   {mean_true_score:.4f}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    precision_recall_path = args.outdir / "precision-recall-plot.png"
    plot_precision_recall(evaluations, precision_recall_path)

    score_correctness_path = args.outdir / "score-correctness-histogram.png"
    plot_score_correctness(evaluations, score_correctness_path)

    print(f"Precision-recall comparison plot saved to {precision_recall_path}")
    print(f"Score correctness histogram saved to {score_correctness_path}")
    


if __name__ == "__main__":
    main()
