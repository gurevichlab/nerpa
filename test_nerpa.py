import argparse
import subprocess
from argparse import Namespace as CommandlineArgs
from itertools import islice
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import polars as pl
import yaml

from src.matching.match_type import Match
from src.pipeline.logging.logger import LoggingConfig, NerpaLogger
from src.testing.check_matches import find_wrong_matches
from src.testing.simplified_alignment import fst_mismatched_step
from src.testing.testing_types import (
    CuratedAlignment,
    TestResult,
    _wrap_alignment,
    simplified_alignment_from_match,
)

NRP_ISO_CLASS = str  # NRP representative ID up to isomorphism
NRP_ID = str
BGC_ID = str

def load_matches(
        nerpa_results_dir: Path,
        nrp_id_to_iso_class: dict[NRP_ID, NRP_ISO_CLASS]
) -> dict[tuple[BGC_ID, NRP_ISO_CLASS], Match]:
    result = {}
    for results_for_bgc in nerpa_results_dir.iterdir():
        if not results_for_bgc.is_dir():
            continue

        with open(results_for_bgc / 'matches_details/matches.yaml') as matches_file:
            loaded_matches: list[Match] = [
                Match.from_dict(match_dict)
                for match_dict in yaml.safe_load(matches_file)
            ]

        for m in loaded_matches:
            bgc_id = m.genome_id
            if bgc_id.endswith('.gbk'):
                bgc_id = bgc_id[:-4]
            iso_class = nrp_id_to_iso_class[m.nrp_variant_id.nrp_id]
            if (bgc_id, iso_class) not in result:
                result[(bgc_id, iso_class)] = m

    return result


def load_command_line_args() -> CommandlineArgs:
    parser = argparse.ArgumentParser(description="Compares the curated alignments with the alignments produced by Nerpa.")
    parser.add_argument("--curated-alignments", type=Path,
                        help=f"Path to the file with curated alignments.")
    parser.add_argument(
        "--nerpa-results-individual",
        type=Path,
        help="The directory with Nerpa results -- each training BGC against its NRPs from the curated alignments."
    )
    parser.add_argument(
        "--nerpa-results-combined",
        type=Path,
        help=(
            "The directory with Nerpa results -- all training BGCs together against MIBiG+Norine. "
            "Used to determine match ranks"
        )
    )
    parser.add_argument(
        "--pnrpdb-info",
        type=Path,
        help=(
            "Path to pnrpdb2_info.tsv file with information about NRPs. "
            "Used to identify isomorphic NRPs."
        )
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
    )

    return parser.parse_args()


def write_wrong_matches(
        wrong_matches: list[dict[str, object]],
        output_file: Path
):
    data = []
    for wrong_match_data in wrong_matches:
        ca: CuratedAlignment = wrong_match_data['curated_alignment']
        nm: Match = wrong_match_data['nerpa_match']
        al_wrong = simplified_alignment_from_match(nm)
        al_right = ca.true_alignment
        mismatched_step_wrong, mismatched_step_right = fst_mismatched_step(al_wrong, al_right)
        mismatched_step_wrong = _wrap_alignment([mismatched_step_wrong]) \
            if mismatched_step_wrong is not None else None
        mismatched_step_right = _wrap_alignment([mismatched_step_right]) \
            if mismatched_step_right is not None else None

        data.append({
            'bgc_id': ca.mibig_id,
            'nrp_id': ca.nrp_id,
            'wrong_alignment': _wrap_alignment(simplified_alignment_from_match(nm)),
            'true_alignment': _wrap_alignment(ca.true_alignment),
            'first_mismatched_step': [
                mismatched_step_wrong,
                mismatched_step_right
            ],
        })
    output_file.write_text(yaml.dump(data, default_flow_style=False, sort_keys=False))


def check_nrp_bgc_variants_parsed(nerpa_results_dir: Path,
                                  bgc_ids: set[str],
                                  nrp_ids: set[str],
                                  nrp_id_to_iso_class: dict[str, str],
                                  ):
    parsed_bgcs_file = nerpa_results_dir / 'preprocessed_input/BGC_variants.yaml'
    parsed_nrps_file = nerpa_results_dir / 'preprocessed_input/NRP_variants.yaml'
    with open(parsed_bgcs_file) as f:
        parsed_bgcs = yaml.safe_load(f)
    parsed_bgc_ids = {Path(bgc_variant['bgc_variant_id']['bgc_id']['antiSMASH_file']).stem
                      for bgc_variant in parsed_bgcs}

    with open(parsed_nrps_file) as f:
        parsed_nrps = yaml.safe_load(f)

    parsed_nrp_ids = {nrp_variant['nrp_variant_id']['nrp_id']
                        for nrp_variant in parsed_nrps}
    parsed_nrp_iso_classes = {nrp_id_to_iso_class[nrp_id]
                              for nrp_id in parsed_nrp_ids}

    missing_bgcs = bgc_ids - parsed_bgc_ids
    missing_nrps = {
        nrp_id
        for nrp_id in nrp_ids
        if nrp_id not in nrp_id_to_iso_class
        or nrp_id_to_iso_class[nrp_id] not in parsed_nrp_iso_classes
    }
    if missing_bgcs:
        # raise ValueError(f'The following BGC IDs are missing in preprocessed input:\n{missing_bgcs}')
        print(f'The following BGC IDs are missing in preprocessed input:\n{missing_bgcs}')

    if missing_nrps:
        # raise ValueError(f'The following NRP IDs are missing (up to isomorphism) in preprocessed input:\n{missing_nrps}')
        print(f'The following NRP IDs are missing (up to isomorphism) in preprocessed input:\n{missing_nrps}')


def remove_deprecated_nrps(test_matches: List[CuratedAlignment],
                           pnrpdb_info: pl.DataFrame) -> List[CuratedAlignment]:
    print('Filtering out matches with deprecated NRPs:')
    '''
    MIBiG is constantly updated and whatnot.
    pnrpdb2_info.tsv is expected to contain the up-to-date list of NRPs.
    I am filtering out test matches with NRP IDs not present in pnrpdb2_info.tsv
    '''

    valid_nrp_ids = set(pnrpdb_info['compound_id'])
    deprecated_ids = {
        test_match.nrp_id
        for test_match in test_matches
        if test_match.nrp_id not in valid_nrp_ids
    }
    print(
        f'Removed {len(deprecated_ids)} test matches with deprecated NRP IDs:\n'
        f'{sorted(deprecated_ids)}\n'
    )

    return [
        test_match
        for test_match in test_matches
        if test_match.nrp_id in valid_nrp_ids
    ]


# TODO: load paths from config instead of hardcoding them


def load_nrp_id_to_iso_class(nerpa_results_dir: Path) -> Dict[NRP_ID, NRP_ISO_CLASS]:
    yaml_file = nerpa_results_dir / 'preprocessed_input/NRP_representatives.yaml'
    with open(yaml_file) as f:
        d = yaml.safe_load(f)

    return {
        entry['nrp_id'][0]: entry['representative_id'][0]
        for entry in d
    }


def load_match_ranks(
        nerpa_results_combined: Path,
        nrp_id_to_iso_class: dict[NRP_ID, NRP_ISO_CLASS]
) -> dict[Tuple[BGC_ID, NRP_ISO_CLASS], int]:
    report = pl.read_csv(nerpa_results_combined / "report.tsv", separator="\t")

    # Rank each NRP against other matches reported for the same BGC.
    ranked_report = (
        report
        .with_columns(
            pl.col("antiSMASH_input")
            .map_elements(lambda path: Path(path).stem, return_dtype=pl.String)
            .alias("BGC_ID")
        )
        .with_columns(
            pl.col("Score")
            .rank(method="ordinal", descending=True)
            .over("BGC_ID")
            .cast(pl.Int64)
            .alias("rank")
        )
    )

    best_ranks = (
        ranked_report
        .with_columns(
            pl.col("NRP_ID")
            .replace_strict(nrp_id_to_iso_class)
            .alias("NRP_ISO_CLASS")
        )
        .group_by("BGC_ID", "NRP_ISO_CLASS")
        .agg(pl.col("rank").min())
    )

    return {
        (bgc_id, iso_class): rank
        for bgc_id, iso_class, rank in best_ranks.iter_rows()
    }


def main():
    nerpa_dir = Path(__file__).parent
    args = load_command_line_args()
    pnrpdb_info: pl.DataFrame = pl.read_csv(args.pnrpdb_info, separator='\t')
    nrp_id_to_iso_class: dict[NRP_ID, NRP_ISO_CLASS] = {
        row['compound_id']: row['nrp_variant_iso_class_representative']
        for row in pnrpdb_info.select(['compound_id', 'nrp_variant_iso_class_representative']).to_dicts()
    }

    # logging_cfg = LoggingConfig(
    #     log_file=args.output_dir / 'test_nerpa.log',
    #     warnings_file=args.output_dir / 'test_nerpa_warnings.log',
    #     stdout_log_level='INFO',
    #     nerpa_version=NERPA_VERSION,
    # )
    #nrp_id_to_iso_class: Dict[NRP_ID, NRP_ISO_CLASS] = {
    #    row['compound_id']: row['nrp_variant_iso_class_representative']
    #    for row in pnrpdb_info.select(['compound_id', 'nrp_variant_iso_class_representative']).to_dicts()
    #}

    # MISSING_EXPECTED: Set[Tuple[BGC_ID, NRP_ID]] = {
    #     (bgc_id, nrp_id)
    #     for bgc_id, nrp_id in
    #     [
    #         ('BGC0001667', 'BGC0001667.4'),   # NRP aligns poorly to its true BGC
    #         ('BGC0000395', 'BGC0000395.0'),  # NRP aligns poorly to its true BGC
    #         ('BGC0001344', 'BGC0001344.8'),  # NRP aligns poorly to its true BGC
    #     ]
    # }


    print('Loading tests')
    curated_alignments_data = yaml.safe_load(args.curated_alignments.read_text(encoding='utf-8'))
    curated_alignments: list[CuratedAlignment] = [
        CuratedAlignment.from_yaml_dict(ca)
        for ca in curated_alignments_data["curated_alignments"]
    ]
    # curated_alignments = remove_deprecated_nrps(curated_alignments, pnrpdb_info)
    curated_alignments_by_id: Dict[Tuple[BGC_ID, NRP_ISO_CLASS], CuratedAlignment] = {
        (test_match.mibig_id, nrp_id_to_iso_class[test_match.nrp_id]): test_match
        for test_match in curated_alignments
    }

    match_ranks: dict[Tuple[BGC_ID, NRP_ISO_CLASS], int] = (
        load_match_ranks(args.nerpa_results_combined, nrp_id_to_iso_class)
    )

    nerpa_matches_by_id: dict[Tuple[BGC_ID, NRP_ISO_CLASS], Match] = (
        load_matches(args.nerpa_results_individual, nrp_id_to_iso_class)
    )

    # check_nrp_bgc_variants_parsed(
    #     nerpa_results_dir=nerpa_results_dir,
    #     bgc_ids={test_match.mibig_id for test_match in curated_alignments},
    #     nrp_ids={test_match.nrp_id for test_match in curated_alignments},
    #     nrp_id_to_iso_class=nrp_id_to_iso_class
    # )
    # matches = load_matches(nerpa_results_dir)
    # matches_by_id: Dict[Tuple[BGC_ID, NRP_ISO_CLASS], Match] = {
    #     (m.genome_id, nrp_id_to_iso_class[m.nrp_variant_id.nrp_id]): m
    #     for m in matches
    # }

    # test_matches_by_id: Dict[Tuple[BGC_ID, NRP_ID], CuratedAlignment] = {
    #     (test_match.mibig_id, test_match.nrp_id): test_match
    #     for test_match in curated_alignments
    # }

    # print('Checking matches')
    # missing_tests = {
    #     (test_match.mibig_id, test_match.nrp_id)
    #     for test_match in curated_alignments
    #     if (test_match.mibig_id, nrp_id_to_iso_class.get(test_match.nrp_id)) not in matches_by_id
    # }

    # missing_expected = missing_tests.intersection(MISSING_EXPECTED)
    # missing_unexpected: Set[Tuple[BGC_ID, NRP_ISO_CLASS]] = missing_tests - missing_expected

    verdicts = {
        (ca.mibig_id, nrp_id_to_iso_class[ca.nrp_id]):
            ca.test(nerpa_matches_by_id[(ca.mibig_id, nrp_id_to_iso_class[ca.nrp_id])],
                    check_ids=False)
        for ca in curated_alignments
        # if (ca.mibig_id, nrp_id_to_iso_class.get(ca.nrp_id)) in matches_by_id
    }

    total = len(curated_alignments)
    verdicts_lst = list(verdicts.values())
    # verdicts_wrong = verdicts_lst.count(TestResult.WRONG)
    # wrong_expected = len({(bgc_id, nrp_id)
    #                       for (bgc_id, nrp_id), verdict in verdicts.items()
    #                       if verdict == TestResult.WRONG and (bgc_id, nrp_id) in WRONG_EXPECTED})
    print(
        f'Correct: {verdicts_lst.count(TestResult.CORRECT)}/{total}\n'
        f'Acceptable alternative: {verdicts_lst.count(TestResult.ACCEPTABLE_ALTERNATIVE)}/{total}\n'
        f'Wrong: {verdicts_lst.count(TestResult.WRONG)}/{total}\n'
    )

    # Missing ranks are treated as falling beyond the reported top 50.
    ranks: list[int] = [
        match_ranks.get(match_id, 51)
        for match_id in verdicts
    ]
    print(
        f"Rank 1: {sum(rank == 1 for rank in ranks)}/{total}\n"
        f"Rank <= 5: {sum(rank <= 5 for rank in ranks)}/{total}\n"
        f"Rank <= 10: {sum(rank <= 10 for rank in ranks)}/{total}\n"
        f"Rank <= 20: {sum(rank <= 20 for rank in ranks)}/{total}\n"
        f"Rank <= 50: {sum(rank <= 50 for rank in ranks)}/{total}\n"
        f"Rank > 50: {sum(rank > 50 for rank in ranks)}/{total}\n"
    )
    

    # print(f'Missing matches (unexpected):\n{sorted(nrp_id for bgc_id, nrp_id in missing_unexpected)}')

    wrong_matches_file = args.output_dir / 'wrong_matches.yaml'
    write_wrong_matches(
        [
            {
                "nerpa_match": nerpa_matches_by_id[(bgc_id, nrp_id_to_iso_class[nrp_id])],
                "curated_alignment": curated_alignments_by_id[(bgc_id, nrp_id_to_iso_class[nrp_id])]
            }
            for (bgc_id, nrp_id), verdict in verdicts.items()
            if verdict == TestResult.WRONG
        ],
        wrong_matches_file
    )


if __name__ == "__main__":
    main()
