from typing import List, Dict, Set
from pathlib import Path
from collections import defaultdict
import argparse
from argparse import Namespace as CommandlineArgs
from logging import Logger

import polars as pl
import yaml

from src.aa_specificity_prediction_model.specificity_prediction_helper import SpecificityPredictionHelper
from src.antismash_parsing.antismash_parser_types import BGC_ID
from src.antismash_parsing.bgc_variant_types import BGC_Variant
from src.build_output.write_results import write_yaml
from src.config import load_monomer_names_helper, load_config, Config
from src.hmm.hmm_constructor.hmm_constructor_state_edge_context_relations import MATCHING_STATE_TYPES
from src.training.logging_config import configure_logging
from src.hmm.hmm_scoring_config import load_hmm_scoring_config
from src.hmm.hmm_scoring_helper import HMMHelper
from src.matching.matcher import NRP_ID
from src.monomer_names_helper import MonomerNamesHelper
from src.rban_parsing.nrp_variant_types import NRP_Variant
from src.testing.testing_types import TestMatch
from src.training.hmm_parameters.extract_data_for_training import extract_data_for_training, DataForTraining
from src.training.hmm_parameters.bgc_alignment_compatibility import bgc_variant_match_compatible
from src.training.hmm_parameters.hmm_infer_edge_params import infer_edge_params
from src.training.hmm_parameters.hmm_infer_emission_params import infer_emission_params
from src.training.hmm_parameters.norine_stats import load_norine_stats
from src.training.hmm_parameters.training_types import MatchWithBGCNRP, EmissionInfo
from src.training.hmm_parameters.write_results import write_params
from test_nerpa import load_local_paths, run_nerpa, remove_deprecated_nrps


def load_all_bgc_variants(bgc_variants_yaml: Path) -> Dict[BGC_ID, List[BGC_Variant]]:  # bgc_id -> bgc_variants
    bgc_variants = defaultdict(list)
    yaml_variants = [BGC_Variant.from_yaml_dict(bgc_variant_dict)
                     for bgc_variant_dict in yaml.safe_load(bgc_variants_yaml.read_text())]
    for bgc_variant in yaml_variants:
        bgc_variants[bgc_variant.bgc_variant_id.bgc_id].append(bgc_variant)
    return bgc_variants


def load_bgc_variants_for_matches(bgc_id_to_bgc_variants: Dict[BGC_ID, List[BGC_Variant]],
                                  matches: List[TestMatch],
                                  log: Logger) -> Dict[NRP_ID, BGC_Variant]:  # nrp_id -> bgc_variant
    nrp_id_to_bgc_variant = {}
    for match in matches:
        nrp_id = match.nrp_id
        genome_id = match.bgc_id
        try:
            # TestMatches come from MIBiG, and math.bgc_id is genome id, e.g. BGC0001234
            # conversely, Nerpa uses tuples (antiSMASH_file, contig_idx, bgc_idx)
            bgc_id = next(bgc_id
                          for bgc_id in bgc_id_to_bgc_variants.keys()
                          if bgc_id.antiSMASH_file.stem == genome_id)

            bgc_variant = next(bgc_variant
                               for bgc_variant in bgc_id_to_bgc_variants[bgc_id]
                               if bgc_variant_match_compatible(bgc_variant, match))
        except StopIteration:
            log.warning(f'Cannot find appropriate BGC variant for {nrp_id}')
            log.debug(f'Match that failed to find variant: {match}')
            continue

        nrp_id_to_bgc_variant[nrp_id] = bgc_variant

    return nrp_id_to_bgc_variant


def dump_emissions_training_data(data_for_training: DataForTraining, output_dir: Path):
    def emission_dict(emission_info: EmissionInfo) -> dict:
        bgc_id, bgc_module, nrp_monomer, state_type = emission_info
        return {'genome': bgc_id.antiSMASH_file,
                'gene': bgc_module.gene_id,
                'a_domain': bgc_module.a_domain_idx,
                'true_residue': nrp_monomer.residue,
                'predictions': {residue: float(score)
                                for residue, score in bgc_module.residue_score.items()},
                'aa34': bgc_module.aa34_code}

    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / 'match_emissions.yaml', 'w') as match_emissions_file:
        yaml.dump([emission_dict(emission_info)
                   for emission_info in data_for_training.emissions
                   if emission_info.state_type in MATCHING_STATE_TYPES],
                  match_emissions_file)


def load_command_line_args(nerpa_dir: Path) -> CommandlineArgs:
    parser = argparse.ArgumentParser(description="Estimates Nerpa parameters "
                                                 "based on Nerpa results on approved matches")
    parser.add_argument("--approved-matches", type=Path,
                        default=nerpa_dir / 'data/for_training_and_testing/approved_matches.yaml')
    parser.add_argument("--bgc-variants", type=Path)
    parser.add_argument("--nrp-variants", type=Path)

    parser.add_argument("--output-dir", type=Path,
                        default=nerpa_dir / 'training_results')
    parser.add_argument("--draw-hmms", action='store_true')
    return parser.parse_args()


def load_hmm_helper(cfg: Config,
                    monomer_names_helper: MonomerNamesHelper) -> HMMHelper:
    specificity_prediction_helper = SpecificityPredictionHelper(cfg.specificity_prediction_config,
                                                                monomer_names_helper)

    hmm_scoring_config = load_hmm_scoring_config(cfg.nerpa_dir,
                                                 cfg.hmm_scoring_config,
                                                 specificity_prediction_helper,
                                                 monomer_names_helper)

    return HMMHelper(hmm_scoring_config, monomer_names_helper)


def get_bgc_variants(bgc_ids: Set[str],
                     nerpa_dir: Path,
                     output_dir: Path) -> Path:
    local_paths = load_local_paths(nerpa_dir)
    antismash_results_dir = local_paths['as_results_mibig4_nrps']

    assert all((antismash_results_dir / bgc_id).exists()
               for bgc_id in bgc_ids), "Some BGC IDs do not have antismash results"

    as_results_file = output_dir / 'antismash_results_paths.txt'
    with open(as_results_file, 'w') as f:
        for bgc_id in bgc_ids:
            f.write(f"{antismash_results_dir / bgc_id}\n")

    smiles_tsv = output_dir / 'smiles.tsv'
    with open(smiles_tsv, 'w') as f:
        f.write('ID\tSMILES\n')
        f.write('compound_1\tCC(=O)OC1=CC=CC=C1C(=O)O\n')  # aspirin as dummy SMILES

    print('Running Nerpa')
    run_nerpa(nerpa_dir,
              antismash_paths=as_results_file,
              smiles_tsv=smiles_tsv,
              min_num_matches_per_bgc=1,
              max_num_matches_per_bgc=1,
              disable_calibration=True,
              disable_dictionary_lookup=True,
              output_dir=output_dir / 'nerpa_results')
    print('Nerpa finished')
    return output_dir / 'nerpa_results/preprocessed_input/BGC_variants.yaml'


def get_nrp_variants(nerpa_dir: Path) -> Dict[NRP_ID, NRP_Variant]:
    nrp_variants_yaml = (nerpa_dir / 'data/input/preprocessed/'
                                     'pnrpdb2_preprocessed.yaml')
    nrp_variants = [NRP_Variant.from_yaml_dict(nrp_variant_dict)
                    for nrp_variant_dict in yaml.safe_load(nrp_variants_yaml.read_text())]
    return {variant.nrp_variant_id.nrp_id: variant
            for variant in nrp_variants}

'''
pnrpdb_info_tsv = nerpa_dir / 'data/for_training_and_testing/pnrpdb2_info.tsv'

pnrpdb_info = pl.read_csv(pnrpdb_info_tsv, separator='\t')

# Create a mapping from nrp_id (compound_id) to nrp_iso_class
nrp_id_to_iso_class = {row['compound_id']: row['nrp_variant_iso_class_representative']
                       for row in pnrpdb_info.iter_rows(named=True)}

# Create a mapping from nrp_iso_class to NRP_Variant
# by finding which iso_class each variant's nrp_id belongs to
iso_class_to_variant = {}
for variant in nrp_variants:
    nrp_id = variant.nrp_variant_id.nrp_id
    iso_class = nrp_id_to_iso_class[nrp_id]
    iso_class_to_variant[iso_class] = variant

# Map each compound_id to its corresponding variant based on nrp_iso_class
nrp_id_to_variant = {}
for row in pnrpdb_info.iter_rows(named=True):
    compound_id = row['compound_id']
    nrp_iso_class = row['nrp_variant_iso_class_representative']

    if nrp_iso_class in iso_class_to_variant:
        nrp_id_to_variant[compound_id] = iso_class_to_variant[nrp_iso_class]
'''

def check_bgcs_with_many_variants(bgc_variants: List[BGC_Variant]):
    bgc_id_to_variants = defaultdict(set)
    for bgc_variant in bgc_variants:
        bgc_id_to_variants[bgc_variant.bgc_variant_id.bgc_id].add(bgc_variant.bgc_variant_id)

    for bgc_id, variant_ids in bgc_id_to_variants.items():
        if len(variant_ids) > 1:
            print(f'Warning: BGC {bgc_id} has {len(variant_ids)} variants used for training:\n'
                  f'{variant_ids}\n')

# TODO: load paths from config instead of hardcoding them
def main():
    nerpa_dir = Path(__file__).parent
    assert nerpa_dir.stem == 'nerpa2', "quick sanity check that nerpa_dir is correct"

    args = load_command_line_args(nerpa_dir)

    logger = configure_logging(args.output_dir)

    # 1. Load configs and helpers
    nerpa_cfg = load_config()
    monomer_names_helper = load_monomer_names_helper(nerpa_cfg.monomers_config,
                                                     nerpa_dir)
    hmm_helper = load_hmm_helper(nerpa_cfg, monomer_names_helper)
    pnrpdb_info = pl.read_csv(nerpa_dir / 'data/for_training_and_testing/pnrpdb2_info.tsv', separator='\t')

    # 2. Load approved matches and corresponding BGC variants
    logger.info('Loading approved matches')
    approved_matches = [TestMatch.from_yaml_dict(test_match_dict)
                        for test_match_dict in yaml.safe_load(args.approved_matches.read_text())]
    approved_matches = remove_deprecated_nrps(approved_matches, pnrpdb_info)

    logger.info('Loading BGC variants')
    # Note: each BGC can have multiple BGC variants;
    # we try to find the ones compatible with alignments in approved matches
    if args.bgc_variants:
        bgc_variants_yaml = args.bgc_variants
    else:
        test_matches_bgcs = {match.bgc_id for match in approved_matches}
        bgc_variants_yaml = get_bgc_variants(test_matches_bgcs,
                                             nerpa_dir,
                                             args.output_dir,)

    bgc_id_to_bgc_variants = load_all_bgc_variants(bgc_variants_yaml)
    matched_bgc_variants: Dict[NRP_ID, BGC_Variant] = load_bgc_variants_for_matches(bgc_id_to_bgc_variants,
                                                                                    approved_matches,
                                                                                    logger)
    check_bgcs_with_many_variants(list(matched_bgc_variants.values()))

    nrp_id_to_nrp_variant = get_nrp_variants(nerpa_dir)


    # note: alignments are taken from "old" approved matches while bgc predictions are taken from current matches
    matches_with_bgcs_nrps_for_training = [
        MatchWithBGCNRP(match,
                        matched_bgc_variants[match.nrp_id],
                        nrp_id_to_nrp_variant[match.nrp_id])
        for match in approved_matches
        if match.nrp_id in matched_bgc_variants
    ]

    logger.info(f'Number of matches with BGC and NRP variants for training: {len(matches_with_bgcs_nrps_for_training)}')

    data_for_training = extract_data_for_training(matches_with_bgcs_nrps_for_training,
                                                  hmm_helper,
                                                  dir_for_hmm_figures=args.output_dir if args.draw_hmms else None)

    # intermediate results for debug
    dump_emissions_training_data(data_for_training, args.output_dir)
    write_yaml(data_for_training.to_dict(), args.output_dir / 'data_for_training.yaml')


    norine_stats = load_norine_stats(nerpa_dir / 'data/nerpa_internal/norine_monomers_info.yaml')
    edge_params = infer_edge_params(data_for_training.edge_choices_per_state_cntxt,
                                    data_for_training.edge_choices_wo_filtering,
                                    data_for_training.edge_choices_df,
                                    logger=logger)

    emission_params = infer_emission_params(data_for_training.emissions,
                                            norine_stats,
                                            args.output_dir,
                                            monomer_names_helper)  # for plotting step function
    write_params(edge_params=edge_params,
                 emission_params=emission_params,
                 data_for_training=data_for_training,
                 old_configs_dir=nerpa_dir / 'configs',
                 output_dir=args.output_dir)


if __name__ == "__main__":
    main()
