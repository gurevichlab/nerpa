from typing import List, Optional, Dict, Tuple
from pathlib import Path
import yaml
from src.matching.match_type import Match
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
from src.matching.detailed_hmm import DetailedHMM
from src.matching.hmm_scoring_helper import HMMHelper
from src.matching.hmm_config import load_hmm_scoring_config
from src.training.paras import substitute_predictions_with_paras
from collections import defaultdict
from src.data_types import BGC_Variant, NRP_Monomer, BGC_Module, Match_BGC_Variant_Info
from src.training.extract_data_for_training import DataForTraining
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
    nerpa_config = yaml.safe_load((nerpa_dir / 'configs/config.yaml').open())
    monomers_cfg = yaml.safe_load((nerpa_dir / nerpa_config['monomers_config']).open('r'))
    monomers_table_tsv = nerpa_dir / monomers_cfg['monomer_names_table']
    return MonomerNamesHelper(pd.read_csv(monomers_table_tsv, sep='\t'),
                              monomers_cfg['supported_residues'],
                              monomers_cfg['pks_names'])

def dump_emissions_training_data(data_for_training: DataForTraining, output_dir: Path):
    def emission_dict(emission_info: Tuple[Match_BGC_Variant_Info, BGC_Module, NRP_Monomer]) -> dict:
        bgc_info, bgc_module, nrp_monomer = emission_info
        return {'genome': bgc_info.genome_id,
                'gene': bgc_module.gene_id,
                'a_domain': bgc_module.a_domain_idx,
                'true_residue': nrp_monomer.residue,
                'predictions': {residue: float(score)
                                for residue, score in bgc_module.residue_score.items()},
                'aa34': bgc_module.aa34_code}

    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / 'match_emissions.yaml', 'w') as match_emissions_file:
        yaml.dump([emission_dict(emission_info)
                   for emission_info in data_for_training.match_emissions],
                  match_emissions_file)


# TODO: load paths from config instead of hardcoding them
def main():
    nerpa_dir = Path(__file__).parent
    monomer_names_helper = load_monomer_names_helper(nerpa_dir)
    hmm_config_file = nerpa_dir / 'configs/hmm_scoring_config.yaml'
    hmm_scoring_config = load_hmm_scoring_config(hmm_config_file)
    hmm_helper = HMMHelper(hmm_scoring_config, monomer_names_helper)

    nerpa_results_dir = nerpa_dir / 'test_results/nerpa_results'

    print('Loading approved matches')
    approved_matches = load_matches_from_txt(nerpa_dir / 'test_data/approved_matches/approved_matches.txt')

    print('Loading BGC variants')
    bgc_variants = load_all_bgc_variants(approved_matches, nerpa_results_dir)

    print('Fixing matches')
    approved_matches = [fix_match(match, monomer_names_helper) for match in approved_matches]

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

    data_for_training = extract_data_for_training(matches_with_bgcs_nrps_for_training,
                                                  hmm_helper)

    # intermediate results for debug
    dump_emissions_training_data(data_for_training, nerpa_dir / 'training_results')

    norine_stats = load_norine_stats(nerpa_dir / 'data/norine_monomers_info.yaml')
    edge_params = infer_edge_params(data_for_training.edge_choices_cnts)
    emission_params = infer_emission_params(data_for_training.match_emissions,
                                            norine_stats,
                                            nerpa_dir / 'training_results')
    write_params(edge_params, emission_params, nerpa_dir / 'training_results')


if __name__ == "__main__":
    main()
