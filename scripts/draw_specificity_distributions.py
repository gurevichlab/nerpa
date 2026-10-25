from __future__ import annotations
import argparse
import json
from math import log
import math
from pathlib import Path
from typing import cast

import matplotlib
import polars as pl
import yaml

from src.config import load_config, load_monomer_names_helper
from src.build_output.draw_graph import RGB, make_color_dict
from src.antismash_parsing.bgc_variant_types import BGC_Module, BGC_Variant
from src.general_type_aliases import Prob


matplotlib.use("Agg")

import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate simplified text and graphical view of specificity distributions for provided BGCs."
    )
    parser.add_argument(
        "--bgc-variants",
        type=Path,
        help="BGC variants yaml file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output directory",
    )
    return parser.parse_args()


def simplify_distribution(
        distribution: dict[str, float],
        max_events: int | None = None,
        min_probability: float = 0.0,
) -> dict[str, float]:
    """
    Leave at most max_events events with probability >= min_probability, and combine the rest into a single "rest" event.
    """
    if max_events is None:
        max_events = len(distribution)
        
    if "rest" in distribution:
        raise ValueError(
            "The input distribution already contains a 'rest' event."
        )
    ranked_events = sorted(
        distribution.items(),
        key=lambda item: item[1],
        reverse=True,
    )
    displayed_events = [
        (event, probability)
        for event, probability in ranked_events
        if probability >= min_probability
    ][:max_events]

    displayed_names = {event for event, _ in displayed_events}
    rest_probability = sum(
        probability
        for event, probability in distribution.items()
        if event not in displayed_names
    )

    simplified_distribution = dict(displayed_events)
    if rest_probability > 0:
        simplified_distribution["rest"] = rest_probability

    return simplified_distribution


class DistributionDisplayHelper:
    distributions: list[dict[str, Prob]]
    names: list[str]
    color_map: dict[str, RGB] | None = None

    def __init__(
            self,
            distributions: list[dict[str, Prob]],
            names: list[str],
            color_map: dict[str, RGB] | None = None,
    ) -> None:
        if len(distributions) != len(names):
            raise ValueError(
                "The number of distributions must match the number of names."
            )
        self.distributions = distributions
        self.names = names
        self.color_map = color_map
        if (color_map is not None
            and 'rest' not in color_map):
            self.color_map['rest'] = (0.5, 0.5, 0.5)  # gray for "rest" event

    @classmethod
    def from_bgc_variant(
            cls,
            bgc_variant: BGC_Variant,
            color_map: dict[str, RGB] | None = None,
    ) -> DistributionDisplayHelper:
        def module_id(module: BGC_Module) -> str:
            return f"{module.gene_id}:{module.a_domain_idx}"

        distributions = []
        names = []
        for module in bgc_variant.modules:
            names.append(module_id(module))
            distributions.append({
                res: math.e ** log_prob
                for res, log_prob in module.residue_score.items()
            })
        return cls(
            distributions=distributions,
            names=names,
            color_map=color_map,
        )

    def to_tsv(
            self,
            output_path: Path,
            max_events: int | None = None,
            min_probability: float = 0.0,
    ) -> None:
        def distribution_to_column(distribution: dict[str, Prob]) -> list[str]:
            d = simplify_distribution(
                distribution,
                max_events=max_events,
                min_probability=min_probability,
            )
            return (
                [
                    f'{res}: {prob:.2f}'
                    for res, prob in sorted(d.items(), key=lambda item: item[1], reverse=True)
                    if res != "rest"
                ]
                + ([f'rest: {d["rest"]:.2f}'] if "rest" in d else [])
            )


        columns = [
            distribution_to_column(distribution)
            for distribution in self.distributions
        ]
        max_len = max(len(column) for column in columns)
        for column in columns:
            column.extend([None] * (max_len - len(column)))
            
        data: dict[str, list[str]] = {
            name: column
            for name, column in zip(self.names, columns, strict=True)
        }

        output_path.parent.mkdir(parents=True, exist_ok=True)
        pl.DataFrame(data).write_csv(output_path, separator="\t")

    def to_svg(
            self,
            output_path: Path,
            max_events: int | None = None,
            min_probability: float = 0.0,
    ) -> None:
            
        figure_width = max(5.0, 4.0 * len(self.distributions))
        figure, axes = plt.subplots(
            1,
            len(self.distributions),
            figsize=(figure_width, 5),
            squeeze=False,
        )

        for axis, name, distribution in zip(
                axes[0],
                self.names,
                self.distributions,
                strict=True,
        ):
            simplified = simplify_distribution(
                distribution,
                max_events=max_events,
                min_probability=min_probability,
            )
            ranked_events = sorted(
                simplified.items(),
                key=lambda item: item[1],
                reverse=True,
            )

            colors = (
                [
                    self.color_map[event]
                    for event, _ in ranked_events
                ]
                if self.color_map is not None
                else None
            )
            probabilities = [
                probability for _, probability in ranked_events
            ]
            _, label_texts = axis.pie(
                probabilities,
                labels=[event for event, _ in ranked_events],
                colors=colors,
                labeldistance=0.6,
                textprops={
                    "fontweight": "bold",
                    "horizontalalignment": "center",
                    "verticalalignment": "center",
                },
                startangle=90,
                counterclock=False,
            )
            for label_text, probability in zip(
                    label_texts,
                    probabilities,
                    strict=True,
            ):
                label_text.set_fontsize(10 + 12 * math.sqrt(probability))
            axis.set_title(name, fontsize=16, fontweight="bold")
            axis.axis("equal")

        output_path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output_path, format="svg", bbox_inches="tight")
        plt.close(figure)

def main() -> None:
    args = parse_args()
            
    nerpa_dir = Path(__file__).parent.parent
    nerpa_cfg = load_config()
    monomer_names_helper = load_monomer_names_helper(
        nerpa_cfg.monomers_config,
        nerpa_dir
    )
    aa_color = make_color_dict(
        labels=monomer_names_helper.supported_residues
    )

    bgc_variants = [
        BGC_Variant.from_dict(bgc_variant_dict)
        for bgc_variant_dict in yaml.safe_load(
                args.bgc_variants.read_text(encoding="utf-8")
        )
    ]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    for bgc_variant in bgc_variants:
        helper = DistributionDisplayHelper.from_bgc_variant(
            bgc_variant,
            color_map=aa_color,
        )
        variant_output_dir: Path = (
            args.output
            / bgc_variant.bgc_variant_id.to_str_short()
        )
        helper.to_tsv(
            output_path=variant_output_dir / "specificity_distributions.tsv",
            max_events=5,
            min_probability=0.05,
        )
        helper.to_svg(
            output_path=variant_output_dir / "specificity_distributions.svg",
            max_events=5,
            min_probability=0.05,
        )


if __name__ == "__main__":
    main()
