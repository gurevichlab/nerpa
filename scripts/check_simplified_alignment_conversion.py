from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import yaml


NERPA_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(NERPA_DIR))

from src.antismash_parsing.bgc_variant_types import BGC_Variant
from src.matching.alignment_type import AlignmentLight
from src.rban_parsing.nrp_variant_types import NRP_Variant
from src.rban_parsing.rban_monomer import rban_name_no_unk_idx
from src.rban_parsing.rban_parser import Parsed_rBAN_Record
from src.rban_parsing.retrieve_nrp_variants import rban_records_to_nrp_variants
from src.testing.simplified_alignment import (
    SimplifiedAlignment,
    SimplifiedAlignmentStep,
    simplified_alignment_to_light_alignments,
    simplified_alignment_to_str,
)
from src.testing.testing_types import CuratedAlignment
from src.training.hmm_parameters.bgc_alignment_compatibility import (
    bgc_variant_match_compatible,
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Check SimplifiedAlignment -> AlignmentLight round trips for all "
            "curated alignments."
        )
    )
    parser.add_argument(
        "--nerpa-results",
        type=Path,
        required=True,
        help=(
            "Nerpa result directory containing "
            "preprocessed_input/BGC_variants.yaml."
        ),
    )
    parser.add_argument(
        "--curated-alignments",
        type=Path,
        default=(
            NERPA_DIR
            / "data"
            / "for_training_and_testing"
            / "curated_alignments.json"
        ),
    )
    parser.add_argument(
        "--parsed-rban-records",
        type=Path,
        default=(
            NERPA_DIR
            / "data"
            / "input"
            / "preprocessed"
            / "pnrpdb2_expanded_parsed_rban_records.yaml"
        ),
    )
    return parser.parse_args()


def load_curated_alignments(path: Path) -> list[CuratedAlignment]:
    """Load true and acceptable alternative alignments from the curated set."""
    data = json.loads(path.read_text(encoding="utf-8"))
    return [
        CuratedAlignment.from_yaml_dict(alignment)
        for alignment in data["curated_alignments"]
    ]


def load_bgc_variants(path: Path) -> dict[str, list[BGC_Variant]]:
    """Load preprocessed BGC variants and group them by MIBiG ID."""
    variants_by_mibig_id: dict[str, list[BGC_Variant]] = defaultdict(list)

    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    for variant_data in data:
        variant = BGC_Variant.from_dict(variant_data)
        mibig_id = variant.bgc_variant_id.bgc_id._genome_id()
        variants_by_mibig_id[mibig_id].append(variant)

    return dict(variants_by_mibig_id)


def load_nrp_variants(path: Path) -> dict[str, NRP_Variant]:
    """Build the NRP variant associated with each preprocessed NRP ID."""
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    records = [
        Parsed_rBAN_Record.from_dict(record_data)
        for record_data in data
    ]

    return {
        variant.nrp_variant_id.nrp_id: variant
        for variant in rban_records_to_nrp_variants(records)
    }


def remove_skips(alignment: SimplifiedAlignment) -> SimplifiedAlignment:
    """Keep matched pairs and normalize unstable unknown-monomer indexes.

    This removes both module skips, whose monomer is ``None``, and NRP
    insertions, whose A-domain identifier is ``None``. Names such as ``X0``
    and ``X1`` are normalized to ``X`` because their indexes are not stable
    identifiers and are ignored by the reconstruction itself.
    """
    return [
        SimplifiedAlignmentStep(
            a_domain_id=step.a_domain_id,
            rban_name=rban_name_no_unk_idx(step.rban_name),
        )
        for step in alignment
        if step.a_domain_id is not None and step.rban_name is not None
    ]


def light_to_simplified(alignment: AlignmentLight) -> SimplifiedAlignment:
    """Discard concrete monomer data while retaining its curated rBAN name."""
    return [
        SimplifiedAlignmentStep(
            a_domain_id=step.a_domain_id,
            rban_name=(
                step.nrp_monomer.rban_name
                if step.nrp_monomer is not None
                else None
            ),
        )
        for step in alignment
    ]


def iter_alignment_cases(
        curated_alignment: CuratedAlignment,
) -> Iterable[tuple[str, SimplifiedAlignment]]:
    """Yield the true alignment followed by every acceptable alternative."""
    yield "true", curated_alignment.true_alignment

    for index, alignment in enumerate(
        curated_alignment.acceptable_alternative_alignments,
        start=1,
    ):
        yield f"alternative {index}", alignment


def main() -> int:
    """Run all reconstruction checks and return a process exit status."""
    arguments = parse_arguments()

    curated_alignments = load_curated_alignments(
        arguments.curated_alignments
    )
    bgc_variants_by_id = load_bgc_variants(
        arguments.nerpa_results
        / "preprocessed_input"
        / "BGC_variants.yaml"
    )
    nrp_variants_by_id = load_nrp_variants(
        arguments.parsed_rban_records
    )

    checked_cases = 0
    checked_reconstructions = 0
    failures = 0

    for curated_alignment in curated_alignments:
        bgc_variant = next(
            (
                variant
                for variant in bgc_variants_by_id.get(
                    curated_alignment.mibig_id,
                    [],
                )
                if bgc_variant_match_compatible(
                    variant,
                    curated_alignment,
                )
            ),
            None,
        )
        nrp_variant = nrp_variants_by_id.get(curated_alignment.nrp_id)

        if bgc_variant is None:
            print(
                f"FAIL {curated_alignment.mibig_id} / "
                f"{curated_alignment.nrp_id}: "
                "no compatible BGC variant"
            )
            failures += 1
            continue

        if nrp_variant is None:
            print(
                f"FAIL {curated_alignment.mibig_id} / "
                f"{curated_alignment.nrp_id}: "
                "NRP variant not found"
            )
            failures += 1
            continue

        for case_name, simplified_alignment in iter_alignment_cases(
            curated_alignment
        ):
            checked_cases += 1
            light_alignments = list(
                simplified_alignment_to_light_alignments(
                    simplified_alignment,
                    bgc_variant,
                    nrp_variant,
                )
            )

            if not light_alignments:
                print(
                    f"FAIL {curated_alignment.mibig_id} / "
                    f"{curated_alignment.nrp_id} / {case_name}: "
                    "no AlignmentLight reconstruction"
                )
                failures += 1
                continue

            expected = remove_skips(simplified_alignment)

            for reconstruction_index, light_alignment in enumerate(
                light_alignments,
                start=1,
            ):
                checked_reconstructions += 1
                reconstructed = remove_skips(
                    light_to_simplified(light_alignment)
                )

                if reconstructed == expected:
                    continue

                failures += 1
                print(
                    f"FAIL {curated_alignment.mibig_id} / "
                    f"{curated_alignment.nrp_id} / {case_name} / "
                    f"reconstruction {reconstruction_index}"
                )
                print("Expected:")
                print(simplified_alignment_to_str(expected), end="")
                print("Reconstructed:")
                print(
                    simplified_alignment_to_str(reconstructed),
                    end="",
                )

    print(
        f"\nChecked {checked_cases} curated alignment cases and "
        f"{checked_reconstructions} reconstructions: "
        f"{failures} failure(s)."
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())