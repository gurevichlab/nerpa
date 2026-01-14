from pathlib import Path
import yaml
from typing import List, Optional, Dict, Tuple, Set
from src.matching.match_type import Match
from src.testing.check_matches import find_wrong_matches
import argparse
from argparse import Namespace as CommandlineArgs
from itertools import islice
import subprocess
import polars as pl
from src.testing.testing_types import (
    TestMatch,
    TestResult,
    simplified_alignment_from_match,
    _wrap_alignment
)
from src.testing.simplified_alignment import fst_mismatched_step


def load_matches(nerpa_results_dir: Path) -> List[Match]:
    with open(nerpa_results_dir / 'matches_details/matches.yaml') as matches_file:
        loaded_matches = [Match.from_dict(match_dict)
                          for match_dict in yaml.safe_load(matches_file)]
    return loaded_matches


def run_nerpa(nerpa_dir: Path,
              antismash_paths: Path,
              smiles_tsv: Path,
              output_dir: Path,
              max_num_matches: int = 0,
              max_num_matches_per_bgc: int = 10,
              min_num_matches_per_bgc: int = 10,
              disable_calibration: bool = False,
              disable_dictionary_lookup: bool = False,
              disable_bgc_deduplication: bool = False,
              disable_nrp_deduplication: bool = False,
              ):
    nerpa_script = nerpa_dir / "nerpa.py"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Construct the command
    command = [
        "python3", str(nerpa_script),
        "--antismash-paths-file", str(antismash_paths),
        "--smiles-tsv", str(smiles_tsv),
        "--col-id", "ID",
        "--output-dir", str(output_dir),
        "--force-output-dir",
        "--max-num-matches", str(max_num_matches),
        "--max-num-matches-per-bgc", str(max_num_matches_per_bgc),
        "--min-num-matches-per-bgc", str(min_num_matches_per_bgc),
        "--skip-molecule-drawing",
        "--threads", "8",
        "--fast-matching",
        "--dump-all-preprocessed",
        "--let-it-crash",
    ]
    if disable_bgc_deduplication:
        command.append("--disable-bgc-deduplication")
    if disable_nrp_deduplication:
        command.append("--disable-nrp-deduplication")
    if disable_calibration:
        command.append("--disable-calibration")
    if disable_dictionary_lookup:
        command.append("--disable-dictionary-lookup")

    # Execute the command and capture output
    print('Executing command:', ' '.join(command))
    try:
        result = subprocess.run(command, check=True)
        print(f"Nerpa run completed successfully")
    except subprocess.CalledProcessError as e:
        print(f"Error running Nerpa:\n{e.stderr}")
        raise


def load_command_line_args(nerpa_dir: Path,
                           local_paths: dict) -> CommandlineArgs:
    parser = argparse.ArgumentParser(description="Runs Nerpa on NRPs and BGCs "
                                                 "from approved matches and checks the results.")
    parser.add_argument("--test-matches", type=Path,
                        default=nerpa_dir / 'data/for_training_and_testing/approved_matches.yaml',
                        help="Path to the approved matches file.")
    parser.add_argument("--smiles-tsv", type=Path,
                        default=nerpa_dir / 'data/input/pnrpdb2_mibig_norine.tsv')
    parser.add_argument("--antismash-results", type=Path,
                        default=local_paths['as_results_mibig4_nrps'])
    parser.add_argument("--output-dir", type=Path,
                        default=nerpa_dir / 'test_results')
    parser.add_argument("--nerpa-results", type=Path)

    return parser.parse_args()


def load_local_paths(nerpa_dir: Path) -> dict:
    with open(nerpa_dir / 'local_paths.yaml') as local_paths_yaml:
        local_paths_dict = yaml.safe_load(local_paths_yaml)
        if any(path is None for path in local_paths_dict.values()):
            raise ValueError(f"Some local paths are not set in {nerpa_dir / 'local_paths.yaml'}"
                             f"Testing cannot be performed.")
        return {key: Path(value) for key, value in local_paths_dict.items()}


def write_wrong_matches(wrong_matches: List[tuple[Match, TestMatch]],
                        output_file: Path):
    data = []
    for match, test_match in wrong_matches:
        al_wrong = simplified_alignment_from_match(match)
        al_right = test_match.true_alignment
        mismatched_step_wrong, mismatched_step_right = fst_mismatched_step(al_wrong, al_right)
        mismatched_step_wrong = _wrap_alignment([mismatched_step_wrong]) \
            if mismatched_step_wrong is not None else None
        mismatched_step_right = _wrap_alignment([mismatched_step_right]) \
            if mismatched_step_right is not None else None

        data.append({
            'bgc_id': test_match.bgc_id,
            'nrp_id': test_match.nrp_id,
            'wrong_alignment': _wrap_alignment(simplified_alignment_from_match(match)),
            'true_alignment': _wrap_alignment(test_match.true_alignment),
            'first_mismatched_step': [mismatched_step_wrong,
                                      mismatched_step_right],
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
    missing_nrps = {nrp_id for nrp_id in nrp_ids
                    if nrp_id_to_iso_class[nrp_id] not in parsed_nrp_iso_classes}
    if missing_bgcs:
        raise ValueError(f'The following BGC IDs are missing in preprocessed input:\n{missing_bgcs}')
    if missing_nrps:
        raise ValueError(f'The following NRP IDs are missing (up to isomorphism) in preprocessed input:\n{missing_nrps}')


def remove_deprecated_nrps(test_matches: List[TestMatch],
                           pnrpdb_info: pl.DataFrame) -> List[TestMatch]:
    print('Filtering out matches with deprecated NRPs:')
    '''
    MIBiG is constantly updated and whatnot.
    pnrpdb2_info.tsv is expected to contain the up-to-date list of NRPs.
    I am filtering out test matches with NRP IDs not present in pnrpdb2_info.tsv
    '''

    valid_nrp_ids = set(pnrpdb_info['compound_id'])
    deprecated_ids = {test_match.nrp_id for test_match in test_matches
                      if test_match.nrp_id not in valid_nrp_ids}
    print(f'Removed {len(deprecated_ids)} test matches with deprecated NRP IDs:\n'
          f'{sorted(deprecated_ids)}\n')

    return [test_match for test_match in test_matches
            if test_match.nrp_id in valid_nrp_ids]

NRP_ISO_CLASS = str  # NRP representative ID up to isomorphism
NRP_ID = str
BGC_ID = str

# TODO: load paths from config instead of hardcoding them


def load_nrp_id_to_iso_class(nerpa_results_dir: Path) -> Dict[NRP_ID, NRP_ISO_CLASS]:
    yaml_file = nerpa_results_dir / 'preprocessed_input/NRP_representatives.yaml'
    with open(yaml_file) as f:
        d = yaml.safe_load(f)

    return {
        entry['nrp_id'][0]: entry['representative_id'][0]
        for entry in d
    }

def main():
    nerpa_dir = Path(__file__).parent
    local_paths = load_local_paths(nerpa_dir)
    pnrpdb_info = pl.read_csv(nerpa_dir / 'data/for_training_and_testing/pnrpdb2_info.tsv',
                              separator='\t')
    #nrp_id_to_iso_class: Dict[NRP_ID, NRP_ISO_CLASS] = {
    #    row['compound_id']: row['nrp_variant_iso_class_representative']
    #    for row in pnrpdb_info.select(['compound_id', 'nrp_variant_iso_class_representative']).to_dicts()
    #}

    MIBIG_ERRORS: Set[BGC_ID] = {
        'BGC0002379',  # the compound is wrong.
    }
    WRONG_EXPECTED: Set[Tuple[BGC_ID, NRP_ID]] = {
        (bgc_id, nrp_id)
        for bgc_id, nrp_id in
        [
            ('BGC0000342', 'BGC0000342.0'),  # Nerpa incorrectly chooses to skip two monomers instead of iterating a gene
        ]
    }
    MISSING_EXPECTED: Set[Tuple[BGC_ID, NRP_ID]] = {
        (bgc_id, nrp_id)
        for bgc_id, nrp_id in
        [
            ('BGC0001667', 'BGC0001667.4'),   # NRP aligns poorly to its true BGC
            ('BGC0000395', 'BGC0000395.0'),  # NRP aligns poorly to its true BGC
            ('BGC0001344', 'BGC0001344.8'),  # NRP aligns poorly to its true BGC
        ]
    }

    args = load_command_line_args(nerpa_dir, local_paths)

    print('Loading tests')
    test_matches_yaml = yaml.safe_load(args.test_matches.read_text())
    test_matches = [TestMatch.from_yaml_dict(test_match_dict)
                    for test_match_dict in test_matches_yaml]
    test_matches = [test_match for test_match in test_matches
                    if test_match.bgc_id not in MIBIG_ERRORS]
    test_matches = remove_deprecated_nrps(test_matches, pnrpdb_info)

    if args.nerpa_results is None:
        # create a file with paths to antismash results for run_nerpa
        test_matches_bgcs = {match.bgc_id for match in test_matches}
        as_results_file = nerpa_dir / args.output_dir / 'antismash_results_paths.txt'
        as_results_file.parent.mkdir(parents=True, exist_ok=True)
        with open(as_results_file, 'w') as f:
            for bgc_id in test_matches_bgcs:
                f.write(f"{args.antismash_results / bgc_id}\n")

        print('Running Nerpa')
        run_nerpa(nerpa_dir,
                  antismash_paths=as_results_file,
                  smiles_tsv=args.smiles_tsv,
                  output_dir=args.output_dir / 'nerpa_results',
                  max_num_matches_per_bgc=30,
                  disable_bgc_deduplication=True,
                  disable_nrp_deduplication=True,)
        print('Nerpa finished')
        nerpa_results_dir = args.output_dir / 'nerpa_results'
    else:
        print(f'Using existing Nerpa results from {args.nerpa_results}')
        nerpa_results_dir = args.nerpa_results
        nerpa_log = (nerpa_results_dir / 'nerpa.log').read_text()
        run_command = next(line for line in nerpa_log.splitlines()
                           if line.startswith('Started with command:')).split('Started with command: ')[1]
        print(f'Nerpa was run with command:\n{run_command}')

    print('Loading Nerpa results')
    nrp_id_to_iso_class = load_nrp_id_to_iso_class(nerpa_results_dir)
    # TODO: remove this line when NRP isomorphism is fully supported in tests
    #  (need to make sure that equal NRP_Monomer with different rBAN names are treated as isomorphic)
    nrp_id_to_iso_class = {nrp_id: nrp_id for nrp_id in nrp_id_to_iso_class}

    check_nrp_bgc_variants_parsed(nerpa_results_dir=nerpa_results_dir,
                                  bgc_ids={test_match.bgc_id for test_match in test_matches},
                                  nrp_ids={test_match.nrp_id for test_match in test_matches},
                                  nrp_id_to_iso_class=nrp_id_to_iso_class)
    matches = load_matches(nerpa_results_dir)
    matches_by_id: Dict[Tuple[BGC_ID, NRP_ISO_CLASS], Match] = {
        (match.genome_id, nrp_id_to_iso_class[match.nrp_variant_id.nrp_id]): match
        for match in matches
    }

    test_matches_by_id: Dict[Tuple[BGC_ID, NRP_ID], TestMatch] = {
        (test_match.bgc_id, test_match.nrp_id): test_match
        for test_match in test_matches
    }

    print('Checking matches')
    missing_tests = {
        (test_match.bgc_id, test_match.nrp_id)
        for test_match in test_matches
        if (test_match.bgc_id, nrp_id_to_iso_class[test_match.nrp_id]) not in matches_by_id
    }

    missing_expected = missing_tests.intersection(MISSING_EXPECTED)
    missing_unexpected: Set[Tuple[BGC_ID, NRP_ISO_CLASS]] = missing_tests - missing_expected

    verdicts = {
        (test_match.bgc_id, test_match.nrp_id):
            test_match.test(matches_by_id.get((test_match.bgc_id,
                                               nrp_id_to_iso_class[test_match.nrp_id])),
                            check_ids=False)
        for test_match in test_matches
        if (test_match.bgc_id, nrp_id_to_iso_class[test_match.nrp_id]) in matches_by_id
    }

    total = len(test_matches)
    verdicts_lst = list(verdicts.values())
    verdicts_wrong = verdicts_lst.count(TestResult.WRONG)
    wrong_expected = len({(bgc_id, nrp_id)
                          for (bgc_id, nrp_id), verdict in verdicts.items()
                          if verdict == TestResult.WRONG and (bgc_id, nrp_id) in WRONG_EXPECTED})
    print(f'Correct: {verdicts_lst.count(TestResult.CORRECT)}/{total}\n'
          f'Acceptable alternative: {verdicts_lst.count(TestResult.ACCEPTABLE_ALTERNATIVE)}/{total}\n'
          f'Wrong (expected): {wrong_expected}/{total}\n'
          f'Wrong (unexpected): {verdicts_wrong - wrong_expected}/{total}\n'
          f'Missing (expected): {len(missing_expected)}/{total}\n'
          f'Missing (unexpected): {len(missing_unexpected)}/{total}\n')

    print(f'Missing matches (unexpected):\n{sorted(nrp_id for bgc_id, nrp_id in missing_unexpected)}')

    wrong_matches_file = args.output_dir / 'wrong_matches.yaml'
    write_wrong_matches(
        [
            (matches_by_id[(bgc_id, nrp_id_to_iso_class[nrp_id])],
             test_matches_by_id[(bgc_id, nrp_id)])
            for (bgc_id, nrp_id), verdict in verdicts.items()
            if verdict == TestResult.WRONG
        ],
        wrong_matches_file
    )


if __name__ == "__main__":
    main()
