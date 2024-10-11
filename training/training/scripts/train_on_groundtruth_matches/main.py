import argparse
import pandas as pd
from pathlib import Path
import yaml
from typing import Dict, List
from run_nerpa import run_nerpa_on_all_pairs
from check_matches import find_wrong_match
from calculate_parameters import calculate_training_parameters
from collections import defaultdict
from write_results import write_results, show_match


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
    parser.add_argument('--norine',
                        type=Path, help='Path to the precomputed norine statistics')  # TODO: accept list of norine monomers and compute stats here
    parser.add_argument('-o', '--output_dir',
                        type=Path, default=Path(__file__).parent / 'output',
                        help='Output directory')
    parser.add_argument('-i', '--ignore_precomputed',
                        action='store_true', help='Ignore precomputed Nerpa results')

    args = parser.parse_args()
    if args.ignore_precomputed:  # for faster debugging
        args.precomputed_nerpa_results = None

    if (args.precomputed_nerpa_results is not None
            and (args.antismash_results_dir is not None or args.rban_results_dir is not None)):
        print('Precomputed Nerpa results provided. antismash and rban results will be ignored')

    return args


def load_matches(nerpa_results_dir: Path) -> List[dict]:
    matches = []
    for results_dir in nerpa_results_dir.iterdir():
        with open(results_dir / 'matches_details/matches.yaml') as matches_file:
            best_match = max(yaml.safe_load(matches_file),
                             key=lambda match: match['NormalisedScore'],
                             default=None)
            if best_match is None:
                with open('warning.txt', 'a') as f:
                    f.write(f'WARNING! {results_dir} has no matches\n')
                continue
            matches.append(best_match)
    return matches


def load_bgc_variants_for_matches(matches: List[dict],
                                  nerpa_results_dir: Path) -> Dict[str, dict]:  # nrp_id -> bgc_variant
    nrp_ids = {match['NRP'] for match in matches}
    bgc_variants = defaultdict(list)
    for nrp_id in nrp_ids:
        bgc_dir = nerpa_results_dir / nrp_id
        try:
            yaml_files = [f for f in (bgc_dir / 'BGC_variants_before_calibration').iterdir()
                          if f.name.endswith('.yaml')]
        except FileNotFoundError:
            print(f'No BGC variants for {nrp_id}')
            raise
        for yaml_file in yaml_files:
            bgc_variants[nrp_id].extend(yaml.safe_load(yaml_file.read_text()))

    nrp_id_to_bgc_variant = {}
    for match in matches:
        try:
            bgc_variant = next(variant for variant in bgc_variants[match['NRP']]
                               if variant['variant_idx'] == match['BGC_variant_idx'])
        except StopIteration:
            print(f'No BGC variant {match["BGC_variant_idx"]} for {match["NRP"]}')
            raise
        nrp_id_to_bgc_variant[match['NRP']] = bgc_variant
    return nrp_id_to_bgc_variant


def load_nrp_variants_for_matches(matches: List[dict],
                                  nerpa_results_dir: Path) -> Dict[str, dict]:  # nrp_id -> nrp_variant
        nrp_variants = {}
        for match in matches:
            nrp_id = match['NRP']
            nrp_id_matches_file = nerpa_results_dir / nrp_id / 'NRP_variants' / f'{nrp_id}.yaml'
            nrp_id_matches = yaml.safe_load(nrp_id_matches_file.read_text())
            nrp_variants[nrp_id] = next(nrp_id_match
                                        for nrp_id_match in nrp_id_matches
                                        if nrp_id_match['variant_idx'] == match['NRP_variant_idx'])
        return nrp_variants


def main():
    args = parse_args()
    matches_table = pd.read_csv(args.matches_table, sep='\t')
    matches_table['NRP variant'] = matches_table['NRP variant'].apply(lambda x: x.split('#')[0])  # TODO: modify the table

    if args.precomputed_nerpa_results is not None:
        nerpa_results_dir = args.precomputed_nerpa_results
    else:
        nerpa_dir = Path('/home/ilianolhin/git/nerpa2')  # TODO: remove hardcoded path
        nerpa_results_dir = run_nerpa_on_all_pairs(matches_table,
                                                   args.antismash_results_dir, args.rban_results_dir,
                                                   nerpa_dir, args.output_dir / 'nerpa_results')

    print('Loading Nerpa results')
    matches = load_matches(nerpa_results_dir)
    bgc_variants = load_bgc_variants_for_matches(matches, nerpa_results_dir)
    nrp_variants = load_nrp_variants_for_matches(matches, nerpa_results_dir)

    nrp_ids_good_matches = list(matches_table[matches_table['Verdict'] == 'good']['NRP variant'])
    approved_matches = [match  # in case some good matches were reconsidered
                        for match in yaml.safe_load(args.approved_matches.read_text())
                        if match['NRP'] in nrp_ids_good_matches]

    print('Checking matches')
    if (wrong_match_info := find_wrong_match(matches, approved_matches)) is not None:
        nrp_id, wrong_match_str = wrong_match_info
        correct_match_str = show_match(next(match for match in approved_matches
                                            if match['NRP'] == nrp_id))
        print(f'Error in matches: {nrp_id}\nAborting')
        (args.output_dir / 'wrong_match.txt').write_text(f'Wrong:\n{wrong_match_str}\n\nCorrect:\n{correct_match_str}')
        print(f'First wrong match is saved in {args.output_dir / "wrong_match.txt"}')
        exit(0)
    else:
        print('All matches are correct')

    matches_with_bgcs_nrps_for_training = [(match, bgc_variants[match['NRP']], nrp_variants[match['NRP']])
                                      for match in matches
                                      if match['NRP'] in nrp_ids_good_matches]
    print('Calculating training parameters')
    # I pass output_dir to save the step function plot. In the future, the function could be made pure
    norine_stats = yaml.safe_load(args.norine.read_text())
    parameters = calculate_training_parameters(matches_with_bgcs_nrps_for_training,
                                               norine_stats, args.output_dir)
    print('Writing results')
    write_results(matches, bgc_variants, nrp_variants,
                  matches_table, parameters, args.output_dir)


if __name__ == "__main__":
    main()