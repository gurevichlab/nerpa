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

def load_matches(nerpa_results_dir: Path) -> List[Match]:
    with open(nerpa_results_dir / 'matches_details/matches.yaml') as matches_file:
        loaded_matches = [Match.from_dict(match_dict)
                          for match_dict in yaml.safe_load(matches_file)]
    return loaded_matches


def load_matches_from_txt(matches_txt: Path) -> List[Match]:
    matches_strs = matches_txt.read_text().split('\n\n')
    matches_strs = [match_str for match_str in matches_strs
                    if match_str.strip()]
    return [Match.from_str(matches_str)
            for matches_str in matches_strs]


def run_nerpa(nerpa_dir: Path,
              antismash_jsons: Path,
              antismash_paths: Path,
              rban_inputs: Path,
              output_dir: Path):
    """
    Run the Nerpa tool with the specified parameters.

    Parameters:
        nerpa_dir (Path): Path to the directory containing nerpa.py.
        antismash_jsons (Path): Path to the antiSMASH input directory.
        rban_inputs (Path): Path to the rBAN input JSON file.
        output_dir (Path): Path to the output directory.
    """
    # Ensure paths are resolved and exist
    nerpa_script = nerpa_dir / "nerpa.py"
    assert nerpa_script.is_file(), f"Nerpa script not found at: {nerpa_script}"
    assert antismash_jsons.is_dir(), f"antiSMASH inputs directory not found: {antismash_jsons}"
    assert rban_inputs.is_file(), f"rBAN inputs JSON not found: {rban_inputs}"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Construct the command
    command = [
        "python3", str(nerpa_script),
        "--antismash", str(antismash_jsons),
        #"--antismash-paths-file", str(antismash_paths),
        "--rban-json", str(rban_inputs),
        #"--paras-results", str(nerpa_dir / "paras/antismash7.1_nrps"),
        "--output-dir", str(output_dir),
        "--force-output-dir",
        "--max-num-matches", "0",
        "--max-num-matches-per-bgc", "10",
        "--debug",
        "--skip-molecule-drawing",
        "--threads", "10",
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


def load_command_line_args(nerpa_dir: Path) -> CommandlineArgs:
    parser = argparse.ArgumentParser(description="Runs Nerpa on NRPs and BGCs "
                                                 "from approved matches and checks the results.")
    parser.add_argument("--approved-matches", type=Path,
                        default=nerpa_dir / 'test_data/approved_matches/approved_matches.txt')
    parser.add_argument("--rban-records", type=Path,
                        default=nerpa_dir / 'test_data/approved_matches/rban_records/merged.json')
    parser.add_argument("--antismash-jsons", type=Path,
                        default=nerpa_dir / 'test_data/approved_matches/antismash_jsons')
    parser.add_argument("--antismash-paths", type=Path,
                        default=nerpa_dir / 'test_data/approved_matches/paths_to_as_uniq.txt')
    parser.add_argument("--output-dir", type=Path,
                        default=nerpa_dir / 'test_results')

    return parser.parse_args()


# TODO: load paths from config instead of hardcoding them
def main():
    nerpa_dir = Path(__file__).parent
    args = load_command_line_args(nerpa_dir)

    print('Loading approved matches')
    approved_matches = load_matches_from_txt(args.approved_matches)

    print('Running Nerpa')
    run_nerpa(nerpa_dir,
              antismash_jsons=args.antismash_jsons,
              antismash_paths=args.antismash_paths,
              rban_inputs=args.rban_records,
              output_dir=args.output_dir / 'nerpa_results')
    print('Nerpa finished')

    print('Loading Nerpa results')
    matches = load_matches(args.output_dir / 'nerpa_results')

    print('Checking matches')
    missing_cnt, wrong_matches = find_wrong_matches(matches, approved_matches)
    print(f'{missing_cnt} matches are missing in Nerpa output')
    print(f'{len(approved_matches) - missing_cnt} matches were tested')

    matches_table_tsv = nerpa_dir / 'matches_inspection' / 'new_matches_inspection_table.tsv'
    matches_table = pd.read_csv(matches_table_tsv, sep='\t')
    matches_table.set_index('NRP variant', inplace=True)

    wrong_matches_should_pass = [(match, correct_match)
                                 for match, correct_match in wrong_matches
                                 if matches_table.at[match.nrp_variant_id.nrp_id, 'Expected Nerpa behaviour']
                                 == 'should pass']
    wrong_matches_can_fail = [match
                              for match, correct_match in wrong_matches
                              if matches_table.at[match.nrp_variant_id.nrp_id, 'Expected Nerpa behaviour']
                              == 'can fail']

    if wrong_matches_can_fail:
        print(f'{len(wrong_matches_can_fail)} tests allowed to fail have failed (thats ok)')
    if wrong_matches_should_pass:
        print(f'{len(wrong_matches_should_pass)} tests that should pass have failed (thats not ok!)')
        wrong_matches_txt = args.output_dir / 'wrong_matches.txt'
        with open(wrong_matches_txt, 'w') as f:
            for nerpa_match, correct_match in wrong_matches_should_pass:
                f.write(f'Nerpa match (wrong):\n{nerpa_match}\n\nCorrect match:\n{correct_match}\n\n')
    else:
        print('All tests that should pass have passed')


if __name__ == "__main__":
    main()
