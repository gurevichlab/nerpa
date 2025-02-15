from typing import List, Optional, Dict
from pathlib import Path
import yaml
from src.matching.matching_types_match import Match
from src.testing.check_matches import find_wrong_matches
from src.data_types import BGC_Variant
from src.training.fix_match import fix_match, fix_bgc_variant, bgc_variant_match_compatible
from src.training.training_types import MatchWithBGCNRP
from src.training.extract_data_for_training import extract_data_for_training
from src.training.hmm_infer_emission_params import infer_emission_params
from src.training.hmm_infer_edge_params import infer_edge_params
from src.training.norine_stats import load_norine_stats
from src.training.write_results import write_params
from src.monomer_names_helper import MonomerNamesHelper
from src.matching.matcher_viterbi_detailed_hmm import DetailedHMM
from src.matching.hmm_scoring_helper import HMMHelper
from src.matching.hmm_config import load_hmm_scoring_config
from collections import defaultdict
import pandas as pd
from itertools import islice
import subprocess


def load_matches_from_txt(matches_txt: Path) -> List[Match]:
    matches_strs = matches_txt.read_text().split('\n\n')
    matches_strs = [match_str for match_str in matches_strs
                    if match_str.strip()]
    return [Match.from_str(matches_str)
            for matches_str in matches_strs]

def load_all_bgc_variants(matches: List[Match],
                          nerpa_results_dir: Path) -> Dict[str, List[BGC_Variant]]:  # bgc_id -> bgc_variants
    nrp_ids = {match.nrp_variant_info.nrp_id for match in matches}
    bgc_variants = defaultdict(list)
    yaml_files = [f for f in (nerpa_results_dir / 'BGC_variants_before_calibration').iterdir()
                  if f.name.endswith('.yaml')]
    for yaml_file in yaml_files:
        yaml_variants = [BGC_Variant.from_yaml_dict(bgc_variant_dict)
                         for bgc_variant_dict in yaml.safe_load(yaml_file.read_text())]
        bgc_id = yaml_variants[0].genome_id
        bgc_variants[bgc_id].extend(yaml_variants)
    return bgc_variants


def load_bgc_variants_for_matches(matches: List[Match],
                                  bgc_variants: Dict[str, List[BGC_Variant]],
                                  check_by_bgc_variant_id: bool = True) -> Dict[str, BGC_Variant]:  # nrp_id -> bgc_variant
    nrp_id_to_bgc_variant = {}
    for match in matches:
        nrp_id = match.nrp_variant_info.nrp_id
        bgc_id = nrp_id.split('.')[0]
        bgc_variant_idx = match.bgc_variant_info.variant_idx
        try:
            if check_by_bgc_variant_id:
                bgc_variant = next(bgc_variant
                                   for bgc_variant in bgc_variants[bgc_id]
                                   if bgc_variant.variant_idx == bgc_variant_idx)
            else:
                bgc_variant = next(bgc_variant
                                   for bgc_variant in bgc_variants[bgc_id]
                                   if bgc_variant_match_compatible(bgc_variant, match))
        except StopIteration:
            print(f'No compatible BGC variant for {nrp_id}')
            with open('warning.txt', 'w') as f:
                f.write(f'WARNING! No compatible BGC variant for {nrp_id}\n')
                f.write(f'Match:\n {match}\n')
            raise
        nrp_id_to_bgc_variant[nrp_id] = bgc_variant
    return nrp_id_to_bgc_variant

def load_monomer_names_helper(nerpa_dir: Path) -> MonomerNamesHelper:
    monomer_names_csv = nerpa_dir / 'configs/monomer_names_table.tsv'
    return MonomerNamesHelper(pd.read_csv(monomer_names_csv, delimiter='\t'))


# TODO: load paths from config instead of hardcoding them
def main():
    nerpa_dir = Path(__file__).parent
    monomer_names_helper = load_monomer_names_helper(nerpa_dir)
    hmm_config_file = nerpa_dir / 'configs/hmm_scoring_config.yaml'
    hmm_scoring_config = load_hmm_scoring_config(hmm_config_file)
    DetailedHMM.hmm_helper = HMMHelper(hmm_scoring_config, monomer_names_helper)

    nerpa_results_dir = nerpa_dir / 'test_results/nerpa_results'

    print('Loading approved matches')
    approved_matches = load_matches_from_txt(nerpa_dir / 'test_data/approved_matches/approved_matches.txt')

    print('Loading BGC variants')
    bgc_variants = load_all_bgc_variants(approved_matches, nerpa_results_dir)

    print('Fixing matches')
    approved_matches = [fix_match(match) for match in approved_matches]

    # As variants ids may differ, load bgc variants based on compatibility, not by id
    bgc_variants_approved_matches = load_bgc_variants_for_matches(approved_matches,
                                                                  bgc_variants,
                                                                  check_by_bgc_variant_id=False)
    bgc_variants_approved_matches = {nrp_id: fix_bgc_variant(bgc_variant)
                                     for nrp_id, bgc_variant in bgc_variants_approved_matches.items()}


    # note: alignments are taken from "old" approved matches while bgc predictions are taken from current matches
    matches_with_bgcs_nrps_for_training = [
        MatchWithBGCNRP(match,
                        bgc_variants_approved_matches[match.nrp_variant_info.nrp_id],
                        None)
        for match in approved_matches
        if match.nrp_variant_info.nrp_id in bgc_variants_approved_matches
    ]

    data_for_training = extract_data_for_training(matches_with_bgcs_nrps_for_training)
    norine_stats = load_norine_stats(nerpa_dir / 'data/norine_monomers_info.yaml')
    edge_params = infer_edge_params(data_for_training.edge_choices_cnts)
    emission_params = infer_emission_params(data_for_training.match_emissions,
                                            norine_stats,
                                            nerpa_dir / 'training_results')
    write_params(edge_params, emission_params, nerpa_dir / 'training_results')


if __name__ == "__main__":
    main()
