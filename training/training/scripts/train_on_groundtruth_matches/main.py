import argparse
import pandas as pd
from pathlib import Path
import yaml
from typing import Dict, List, NamedTuple, Optional
from run_nerpa import run_nerpa_on_all_pairs
from check_matches import find_wrong_match
#from calculate_parameters import calculate_training_parameters
from collections import defaultdict
from src.matching.matching_types_match import Match
from src.data_types import BGC_Variant, NRP_Variant
from src.write_results import write_yaml
from src.monomer_names_helper import MonomerNamesHelper
from src.monomer_names_helper import monomer_names_helper as mon_helper
from src.matching.matcher_viterbi_detailed_hmm import DetailedHMM
from src.matching.hmm_config import load_hmm_scoring_config
from src.matching.hmm_scoring_helper import HMMHelper
from fix_match import bgc_variant_match_compatible, fix_match, fix_bgc_variant
from norine_stats import NorineStats
from auxilary_types import MatchWithBGCNRP
from itertools import islice, chain
import dacite
from train_hmm import estimate_parameters


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


def load_matches(nerpa_results_dir: Path,
                 max_num_matches: Optional[int] = None) -> List[Match]:
    matches = []
    for results_dir in islice(nerpa_results_dir.iterdir(), max_num_matches):
        with open(results_dir / 'matches_details/matches.yaml') as matches_file:
            loaded_matches = map(lambda match_dict: Match.from_dict(match_dict),
                                 yaml.safe_load(matches_file))
            best_match = max(loaded_matches,
                             key=lambda match: match.normalized_score,
                             default=None)
            if best_match is None:
                with open('warning.txt', 'a') as f:
                    f.write(f'WARNING! {results_dir} has no matches\n')
                continue
            matches.append(best_match)
    return matches


def load_bgc_variants_for_matches(matches: List[Match],
                                  nerpa_results_dir: Path,
                                  check_by_bgc_variant_id: bool = True) -> Dict[str, BGC_Variant]:  # nrp_id -> bgc_variant
    nrp_ids = {match.nrp_variant_info.nrp_id for match in matches}
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
            bgc_variants[nrp_id].extend(BGC_Variant.from_yaml_dict(bgc_variant_dict)
                                        for bgc_variant_dict in yaml.safe_load(yaml_file.read_text()))

    nrp_id_to_bgc_variant = {}
    for match in matches:
        nrp_id = match.nrp_variant_info.nrp_id
        bgc_variant_idx = match.bgc_variant_info.variant_idx
        try:
            if check_by_bgc_variant_id:
                bgc_variant = next(bgc_variant
                                   for bgc_variant in bgc_variants[nrp_id]
                                   if bgc_variant.variant_idx == bgc_variant_idx)
            else:
                bgc_variant = next(bgc_variant
                                   for bgc_variant in bgc_variants[nrp_id]
                                   if bgc_variant_match_compatible(bgc_variant, match))
        except StopIteration:
            print(f'No BGC variant {bgc_variant_idx} for {nrp_id}')
            continue
        nrp_id_to_bgc_variant[nrp_id] = bgc_variant
    return nrp_id_to_bgc_variant


def load_nrp_variants_for_matches(matches: List[Match],
                                  nerpa_results_dir: Path) -> Dict[str, NRP_Variant]:  # nrp_id -> nrp_variant
        nrp_variants = {}
        for match in matches:
            nrp_id = match.nrp_variant_info.nrp_id
            yaml_file = nerpa_results_dir / nrp_id / 'NRP_variants' / f'{nrp_id}.yaml'
            loaded_nrp_variants = (NRP_Variant.from_yaml_dict(nrp_variant_dict)
                                   for nrp_variant_dict in yaml.safe_load(yaml_file.read_text()))
            nrp_variants[nrp_id] = next(nrp_variant
                                        for nrp_variant in loaded_nrp_variants
                                        if nrp_variant.variant_idx == match.nrp_variant_info.variant_idx)
        return nrp_variants


def load_matches_from_txt(matches_txt: Path) -> List[Match]:
    matches_strs = matches_txt.read_text().split('\n\n')
    matches_strs = [match_str for match_str in matches_strs
                    if match_str.strip()]
    return [Match.from_str(matches_str)
            for matches_str in matches_strs]

def load_new_approved_matches(nrp_ids_good_matches: List[str]) -> List[Match]:
    matches_txt = Path('/home/ilianolhin/git/nerpa2/training/training/matches_for_inspection/matches_to_investigate_many_fragments.txt')
    return [match for match in load_matches_from_txt(matches_txt)
            if match.nrp_variant_info.nrp_id in nrp_ids_good_matches]


def load_corrected_approved_matches(nrp_ids_good_matches: List[str]) -> List[Match]:
    matches_txt = Path('/home/ilianolhin/git/nerpa2/training/training/scripts/train_on_groundtruth_matches/corrected_matches')
    return [match for match in load_matches_from_txt(matches_txt)
            if match.nrp_variant_info.nrp_id in nrp_ids_good_matches]


def load_approved_matches(approved_matches: Path, nrp_ids_good_matches: List[str]) -> List[Match]:

    old_approved_matches = [Match.from_dict(match_dict)
                            for match_dict in yaml.safe_load(approved_matches.read_text())]
    corrected_approved_matches = load_corrected_approved_matches(nrp_ids_good_matches)
    new_approved_matches = load_new_approved_matches(nrp_ids_good_matches)

    approved_matches = []
    seen_nrp_ids = set()
    for match in chain(corrected_approved_matches, new_approved_matches, old_approved_matches):
        nrp_id = match.nrp_variant_info.nrp_id
        if nrp_id in nrp_ids_good_matches and nrp_id not in seen_nrp_ids:
            approved_matches.append(match)
            seen_nrp_ids.add(match.nrp_variant_info.nrp_id)
    return approved_matches


def load_monomer_names_helper() -> MonomerNamesHelper:
    monomer_names_csv = Path('/home/ilianolhin/git/nerpa2/configs/monomer_names_table.tsv')
    return MonomerNamesHelper(pd.read_csv(monomer_names_csv, delimiter='\t'))

def main():
    max_num_matches = None  # for debugging
    monomer_names_helper = load_monomer_names_helper()
    hmm_scoring_config = load_hmm_scoring_config(Path('/home/ilianolhin/git/nerpa2/configs/hmm_scoring_config.yaml'))
    DetailedHMM.hmm_helper = HMMHelper(hmm_scoring_config, monomer_names_helper)

    args = parse_args()
    matches_table = pd.read_csv(args.matches_table, sep='\t')
    matches_table['NRP variant'] = matches_table['NRP variant'].apply(lambda x: x.split('#')[0])  # TODO: modify the table
    nrp_ids_good_matches = list(matches_table[matches_table['Verdict'] == 'good']['NRP variant'])
    #nrp_ids_good_matches = ['BGC0001620.1']
    if args.precomputed_nerpa_results is not None:
        nerpa_results_dir = args.precomputed_nerpa_results
    else:
        nerpa_dir = Path('/home/ilianolhin/git/nerpa2')  # TODO: remove hardcoded path
        nerpa_results_dir = run_nerpa_on_all_pairs(matches_table,
                                                   args.antismash_results_dir, args.rban_results_dir,
                                                   nerpa_dir, args.output_dir / 'nerpa_results',
                                                   max_num_pairs=max_num_matches)

    # currently it is a mess but in the future it will be just one file
    print('Loading approved matches')
    approved_matches = load_approved_matches(args.approved_matches, nrp_ids_good_matches)
    print('Fixing matches')
    approved_matches = [fix_match(match) for match in approved_matches]
    # As variants ids may differ, load bgc variants based on compatibility, not by id
    bgc_variants_approved_matches = load_bgc_variants_for_matches(approved_matches,
                                                                  nerpa_results_dir,
                                                                  check_by_bgc_variant_id=False)
    bgc_variants_approved_matches = {nrp_id: fix_bgc_variant(bgc_variant)
                                     for nrp_id, bgc_variant in bgc_variants_approved_matches.items()}

    '''
    print('Loading Nerpa results')
    matches = load_matches(nerpa_results_dir, max_num_matches=max_num_matches)
    bgc_variants = load_bgc_variants_for_matches(matches, nerpa_results_dir)
    nrp_variants = load_nrp_variants_for_matches(matches, nerpa_results_dir)


    print('Checking matches')
    check = False
    if check:
        if (wrong_match_info := find_wrong_match(matches, approved_matches)) is not None:
            wrong_match, correct_match = wrong_match_info
            print(f'Error in matches: {wrong_match.nrp_variant_info.nrp_id}\nAborting')
            (args.output_dir / 'wrong_match.txt').write_text(f'Wrong:\n{wrong_match}\n\nCorrect:\n{correct_match}')
            print(f'First wrong match is saved in {args.output_dir / "wrong_match.txt"}')
            exit(0)
        else:
            print('All matches are correct')
    '''

    # note: alignments are taken from "old" approved matches while bgc predictions are taken from current matches
    matches_with_bgcs_nrps_for_training = [
        MatchWithBGCNRP(match,
                        bgc_variants_approved_matches[match.nrp_variant_info.nrp_id],
                        None)
        for match in approved_matches
        if match.nrp_variant_info.nrp_id in bgc_variants_approved_matches
    ]
    '''
    matches_with_bgcs_nrps_for_training = [MatchWithBGCNRP(match,
                                                           bgc_variants_approved_matches[match.nrp_variant_info.nrp_id],
                                                           nrp_variants[match.nrp_variant_info.nrp_id])
                                           for match in approved_matches
                                           if match.nrp_variant_info.nrp_id in bgc_variants_approved_matches
                                             and match.nrp_variant_info.nrp_id in nrp_variants
                                           ]
    '''

    with open('matches_for_training.txt', 'w') as out:
        out.write('\n\n'.join(str(match) for match, _, _ in matches_with_bgcs_nrps_for_training))
    print('Calculating training parameters')
    hmm_parameters = estimate_parameters(matches_with_bgcs_nrps_for_training)

    # I need to convert the nested dict to plain dict to be able to save it to yaml
    hmm_params_plain_dict = {}
    for edge_type, context_to_prob in hmm_parameters.items():
        context_to_prob_plain = {context.name if context is not None else 'None': prob for context, prob in context_to_prob.items()}
        hmm_params_plain_dict[edge_type.name] = context_to_prob_plain

    write_yaml(hmm_params_plain_dict, args.output_dir / 'hmm_parameters.yaml')
    # I pass output_dir to save the step function plot. In the future, the function could be made pure
    #norine_stats = dacite.from_dict(NorineStats, yaml.safe_load(args.norine.read_text()))
    #parameters = calculate_training_parameters(matches_with_bgcs_nrps_for_training,
    #                                           norine_stats, args.output_dir)
    print('Writing results')
    #write_results(matches, bgc_variants, nrp_variants,
    #              matches_table, parameters, args.output_dir)


if __name__ == "__main__":
    main()