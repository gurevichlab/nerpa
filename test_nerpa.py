from pathlib import Path
import yaml
from typing import List, Optional
from src.matching.match_type import Match
from src.testing.check_matches import find_wrong_matches
import argparse
from argparse import Namespace as CommandlineArgs
from itertools import islice
import subprocess
import pandas as pd
from src.testing.testing_types import (
    TestMatch,
    TestResult,
    simplified_alignment_from_match,
    _wrap_alignment,
)


def load_matches(nerpa_results_dir: Path) -> List[Match]:
    with open(nerpa_results_dir / 'matches_details/matches.yaml') as matches_file:
        loaded_matches = [Match.from_dict(match_dict)
                          for match_dict in yaml.safe_load(matches_file)]
    return loaded_matches


def run_nerpa(nerpa_dir: Path,
              antismash_paths: Path,
              smiles_tsv: Path,
              output_dir: Path):
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
        "--max-num-matches", "0",
        "--max-num-matches-per-bgc", "10",
        "--skip-molecule-drawing",
        "--threads", "8",
        "--fast-matching"
    ]

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
                        default=nerpa_dir / 'data/input/pnrpdb2_mibig_norine_deduplicated.tsv')
    parser.add_argument("--antismash-results", type=Path,
                        default=local_paths['as_results_mibig4_nrps'])
    parser.add_argument("--output-dir", type=Path,
                        default=nerpa_dir / 'test_results')

    return parser.parse_args()


def load_local_paths(nerpa_dir: Path) -> dict:
    with open(nerpa_dir / 'local_paths.yaml') as local_paths_yaml:
        return yaml.safe_load(local_paths_yaml)


def write_wrong_matches(wrong_matches: List[tuple[Match, TestMatch]],
                        output_file: Path):
    data = []
    for match, test_match in wrong_matches:
        data.append({
            'bgc_id': test_match.bgc_id,
            'nrp_id': test_match.nrp_id,
            'wrong_alignment': _wrap_alignment(simplified_alignment_from_match(match)),
            'true_alignment': _wrap_alignment(test_match.true_alignment),
        })
    output_file.write_text(yaml.dump(data, default_flow_style=False, sort_keys=False))


# TODO: load paths from config instead of hardcoding them
def main():
    nerpa_dir = Path(__file__).parent
    local_paths = load_local_paths(nerpa_dir)
    args = load_command_line_args(nerpa_dir, local_paths)

    print('Loading tests')
    test_matches_yaml = yaml.safe_load(args.test_matches.read_text())
    test_matches = [TestMatch.from_yaml_dict(test_match_dict)
                    for test_match_dict in test_matches_yaml]
    test_matches_bgcs = {match.bgc_id for match in test_matches}

    # create a file with paths to antismash results for run_nerpa
    as_results_file = nerpa_dir / 'antismash_paths_for_testing.txt'
    with open(as_results_file, 'w') as f:
        for bgc_id in test_matches_bgcs:
            f.write(f"{args.antismash_results / bgc_id}\n")

    print('Running Nerpa')
    '''
    run_nerpa(nerpa_dir,
              antismash_paths=as_results_file,
              smiles_tsv=args.smiles_tsv,
              output_dir=args.output_dir / 'nerpa_results')
    '''
    print('Nerpa finished')

    print('Loading Nerpa results')
    matches = load_matches(args.output_dir / 'nerpa_results')
    matches_by_id = {(match.bgc_variant_id.bgc_id.genome_id,
                      match.nrp_variant_id.nrp_id): match
                         for match in matches}
    test_matches_by_id = {(test_match.bgc_id, test_match.nrp_id): test_match
                          for test_match in test_matches}


    print('Checking matches')
    missing_tests = test_matches_by_id.keys() - matches_by_id.keys()

    verdicts = {
        (test_match.bgc_id, test_match.nrp_id):
            test_match.test(matches_by_id.get((test_match.bgc_id,
                                               test_match.nrp_id)))
        for test_match in test_matches
        if (test_match.bgc_id, test_match.nrp_id) in matches_by_id
    }

    total = len(test_matches)
    verdicts_lst = list(verdicts.values())
    print(f'Ok: {verdicts_lst.count(TestResult.CORRECT)}/{total}\n'
          f'Acceptable alternative: {verdicts_lst.count(TestResult.ACCEPTABLE_ALTERNATIVE)}/{total}\n'
          f'Wrong: {verdicts_lst.count(TestResult.WRONG)}/{total}\n'
          f'Missing: {len(missing_tests)}/{total}')

    print(f'Missing matches:\n{missing_tests}')

    wrong_matches_file = args.output_dir / 'wrong_matches.yaml'
    write_wrong_matches([(matches_by_id[(bgc_id, nrp_id)],
                          test_matches_by_id[(bgc_id, nrp_id)])
                         for (bgc_id, nrp_id), verdict in verdicts.items()
                         if verdict == TestResult.WRONG],
                       wrong_matches_file)



if __name__ == "__main__":
    main()
