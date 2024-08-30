import argparse
import pandas as pd
from pathlib import Path
import yaml


def parse_args():
    parser = argparse.ArgumentParser(description='Train NERPA on groundtruth matches')
    parser.add_argument('-t', '--matches_table',
                        type=Path, help='Path to the matches table')
    parser.add_argument('-m', '--approved_matches',
                        type=Path, help='Path to the approved matches yaml')
    parser.add_argument('-n', '--precomputed_nerpa_results',
                        type=Path, help='Path to the precomputed nerpa results')
    parser.add_argument('-a', '--antismash_results_dir',
                        type=Path, help='Path to the antismash results on MIBiG')
    parser.add_argument('-r', '--rban_results_dir',
                        type=Path, help='Path to the rban results on MIBiG')
    parser.add_argument('-o', '--output_dir',
                        type=Path, default=Path(__file__).parent / 'output',
                        help='Output directory')

    args = parser.parse_args()
    if args.n is not None and (args.a is not None or args.r is not None):
        parser.error('Precomputed NERPA results should be provided without antismash and rban results')

    return args


def load_matches(nerpa_results_dir: Path) -> List[dict]:
    matches = []
    for results_dir in nerpa_results_dir.iterdir():
        with open(results_dir / 'matches_details/matches.yaml') as matches_file:
            best_match = max(yaml.safe_load(matches_file),
                             key=lambda match: match['NormalisedScore'],
                             default=None)
            if best_match is None:
                print(f'{results_dir} has no matches')
                raise SystemExit
            matches.append(best_match)
    return matches


def main():
    args = parse_args()
    matches_table = pd.read_csv(args.matches_table, sep='\t')
    matches_table['NRP variant'] = matches_table['NRP variant'].apply(lambda x: x.split('#')[0])  # TODO: modify the table

    if args.precomputed_nerpa_results is not None:
        nerpa_results_dir = args.precomputed_nerpa_results
    else:
        nerpa_results_dir = run_nerpa_on_groundtruth_matches(matches_table, args.a, args.r, args.output_dir)

    matches = load_matches(nerpa_results_dir)
    approved_matches = yaml.safe_load(args.approved_matches.read_text())

    all_correct, error_info = check_matches(matches, approved_matches, args.output_dir)
    if all_correct:
        print('All matches are correct')
    else:
        print(f'Error in matches: {error_info.wrong_field}\nAborting')
        (args.output_dir / 'wrong_match.txt').write_text(error_info.wrong_match)
        print(f'First wrong match is saved in {args.output_dir / "wrong_match.txt"}')
        exit(0)

    parameters = calculate_parameters(matches, approved_matches)
    write_results(matches, approved_matches, parameters, args.output_dir)


if __name__ == "__main__":
    main()