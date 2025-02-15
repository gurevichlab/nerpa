from pathlib import Path
import yaml
from typing import List, Optional
from src.matching.match_type import Match
from src.testing.check_matches import find_wrong_matches
from itertools import islice
import subprocess

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


def run_nerpa(nerpa_dir: Path, antismash_inputs: Path, rban_inputs: Path, output_dir: Path):
    """
    Run the Nerpa tool with the specified parameters.

    Parameters:
        nerpa_dir (Path): Path to the directory containing nerpa.py.
        antismash_inputs (Path): Path to the antiSMASH input directory.
        rban_inputs (Path): Path to the rBAN input JSON file.
        output_dir (Path): Path to the output directory.
    """
    # Ensure paths are resolved and exist
    nerpa_script = nerpa_dir / "nerpa.py"
    assert nerpa_script.is_file(), f"Nerpa script not found at: {nerpa_script}"
    assert antismash_inputs.is_dir(), f"antiSMASH inputs directory not found: {antismash_inputs}"
    assert rban_inputs.is_file(), f"rBAN inputs JSON not found: {rban_inputs}"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Construct the command
    command = [
        "python3", str(nerpa_script),
        #"--antismash", str(antismash_inputs),
        "--bgc-variants", "/home/ilianolhin/git/nerpa2/input_all/",
        "--rban-json", str(rban_inputs),
        "--output_dir", str(output_dir),
        "--force-existing-outdir",
        "--max-num-matches", "0",
        "--max-num-matches-per-bgc", "10",
        "--dont-draw-molecules",
        "--threads", "10",
    ]

    # Execute the command and capture output
    print('Executing command:', ' '.join(command))
    try:
        result = subprocess.run(command, check=True)
        print(f"Nerpa run completed successfully. Output:\n{result.stdout}")
    except subprocess.CalledProcessError as e:
        print(f"Error running Nerpa:\n{e.stderr}")
        raise


# TODO: load paths from config instead of hardcoding them
def main():
    nerpa_dir = Path(__file__).parent

    print('Loading approved matches')
    approved_matches = load_matches_from_txt(nerpa_dir / 'test_data/approved_matches/approved_matches.txt')

    print('Running Nerpa')
    run_nerpa(nerpa_dir,
              antismash_inputs=nerpa_dir / 'test_data/approved_matches/antismash_jsons',
              rban_inputs=nerpa_dir / 'test_data/approved_matches/rban_records/merged.json',
              output_dir=nerpa_dir / 'test_results/nerpa_results')
    print('Nerpa finished')

    print('Loading Nerpa results')
    matches = load_matches(nerpa_dir / 'test_results/nerpa_results')

    print('Checking matches')
    check_results = list(find_wrong_matches(matches, approved_matches))

    if check_results:
        print(f'{len(check_results)}/{len(approved_matches)} matches are incorrect')
        wrong_matches_txt = nerpa_dir / 'test_results/wrong_matches.txt'
        with open(wrong_matches_txt, 'w') as f:
            for nerpa_match, correct_match in check_results:
                f.write(f'Nerpa match (wrong):\n{nerpa_match}\n\nCorrect match:\n{correct_match}\n\n')
    else:
        print('All matches are correct')


if __name__ == "__main__":
    main()
