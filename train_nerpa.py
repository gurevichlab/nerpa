import argparse
from argparse import Namespace as CommandlineArgs
from collections import defaultdict
from logging import Logger
from pathlib import Path
from typing import Dict, List, NamedTuple, Set, Optional

import polars as pl
import yaml

from src.aa_specificity_prediction_model.specificity_prediction_helper import (
    SpecificityPredictionHelper,
)
from src.antismash_parsing.antismash_parser_types import BGC_ID
from src.antismash_parsing.bgc_variant_types import BGC_Variant
from src.build_output.write_results import write_yaml
from src.config import Config, load_config, load_monomer_names_helper
from src.hmm.hmm_constructor.hmm_constructor_state_edge_context_relations import (
    MATCHING_STATE_TYPES,
)
from src.hmm.hmm_scoring_config import load_hmm_scoring_config
from src.hmm.hmm_scoring_helper import HMMHelper
from src.matching.matcher import NRP_ID
from src.monomer_names_helper import MonomerNamesHelper
from src.rban_parsing.nrp_variant_types import NRP_Variant
from src.rban_parsing.rban_parser import Parsed_rBAN_Record
from src.rban_parsing.retrieve_nrp_variants import rban_records_to_nrp_variants
from src.testing.testing_types import TestMatch
from src.training.hmm_parameters.bgc_alignment_compatibility import (
    bgc_variant_match_compatible,
)
from src.training.hmm_parameters.extract_data_for_training import (
    DataForTraining,
    extract_data_for_training,
)
from src.training.hmm_parameters.hmm_infer_edge_params import infer_edge_params
from src.training.hmm_parameters.hmm_infer_emission_params import infer_emission_params
from src.training.hmm_parameters.norine_stats import load_norine_stats
from src.training.hmm_parameters.training_types import EmissionInfo, MatchWithBGCNRP
from src.training.hmm_parameters.write_results import write_params
from src.training.logging_config import configure_logging
from test_nerpa import load_local_paths, remove_deprecated_nrps, run_nerpa


def load_all_bgc_variants(bgc_variants_yaml: Path) -> Dict[BGC_ID, List[BGC_Variant]]:  # bgc_id -> bgc_variants
    bgc_variants = defaultdict(list)
    yaml_variants = [BGC_Variant.from_yaml_dict(bgc_variant_dict)
                     for bgc_variant_dict in yaml.safe_load(bgc_variants_yaml.read_text())]
    for bgc_variant in yaml_variants:
        bgc_variants[bgc_variant.bgc_variant_id.bgc_id].append(bgc_variant)
    return bgc_variants


def load_bgc_variants_for_matches(
        bgc_id_to_bgc_variants: Dict[BGC_ID, List[BGC_Variant]],
        matches: List[TestMatch],
        log: Logger
) -> Dict[NRP_ID, BGC_Variant]:  # nrp_id -> bgc_variant
    nrp_id_to_bgc_variant = {}
    for m in matches:
        nrp_id = m.nrp_id
        genome_id = m.bgc_id
        try:
            # TestMatches come from MIBiG, and math.bgc_id is genome id, e.g. BGC0001234
            # conversely, Nerpa uses tuples (antiSMASH_file, contig_idx, bgc_idx)
            bgc_id = next(bgc_id
                          for bgc_id in bgc_id_to_bgc_variants.keys()
                          if bgc_id._genome_id() == genome_id)

            bgc_variant = next(bgc_variant
                               for bgc_variant in bgc_id_to_bgc_variants[bgc_id]
                               if bgc_variant_match_compatible(bgc_variant, m))
        except StopIteration:
            log.warning(f'Cannot find appropriate BGC variant for {nrp_id}')
            log.debug(f'Match that failed to find variant: {m}')
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
    default_approved_matches = (
        nerpa_dir
        / 'data'
        / 'for_training_and_testing'
        / 'approved_matches.yaml'
    )
    parser.add_argument(
        "--approved-matches",
        type=Path,
        default=default_approved_matches,
        help=f"Path to the approved matches YAML file (default: {default_approved_matches})",
    )

    parser.add_argument(
        "--nerpa-results-on-mibig-no-calibration",
        type=Path,
        required=True,
        help=(
            "Path to the directory containing Nerpa results on MIBiG BGCs. "
            "Needed to extract preprocessed BGC variants for training. "
            "PARAS predictions should NOT be calibrated."
        ),
    )

    default_pnrpdb2_info = (
        nerpa_dir
        / 'data'
        / 'for_training_and_testing'
        / 'pnrpdb2_info.tsv'
    )
    parser.add_argument(
        "--pnrpdb2-info",
        type=Path,
        default=default_pnrpdb2_info,
        help=(
            f"Path to the PNRPDB2 info TSV file (default: {default_pnrpdb2_info}). "
            "This file contains compounds clustering and some status such as number of monomers, number of unknown monomers, etc. "
        ),
    )

    default_mibig_norine_preprocessed = (
        nerpa_dir
        / 'data'
        / 'input'
        / 'preprocessed'
        / 'pnrpdb2_mibig_norine_deduplicated_parsed_rban_records.yaml'
    )
    parser.add_argument(
        "--mibig-norine-preprocessed",
        type=Path,
        default=default_mibig_norine_preprocessed,
        help=f"Path to the preprocessed PNRPDB2 data (default: {default_mibig_norine_preprocessed})",
    )
        
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


def get_bgc_variants(
        nerpa_results: Path,
        genome_ids: Optional[Set[str]]=None,
) -> List[BGC_Variant]:
    bgc_variants_yaml = nerpa_results / 'preprocessed_input' / 'BGC_variants.yaml'
    with open(bgc_variants_yaml, 'r') as f:
        bgc_variants = [
            BGC_Variant.from_yaml_dict(bgc_variant_dict)
            for bgc_variant_dict in yaml.safe_load(f)
        ]
    return [
        bgc_variant for bgc_variant in bgc_variants
        if genome_ids is None or bgc_variant.bgc_variant_id.bgc_id._genome_id() in genome_ids
    ]


 def get_nrp_variants(parsed_rban_records_yaml: Path) -> Dict[NRP_ID, NRP_Variant]:
     with open (parsed_rban_records_yaml, 'r') as f:
         parsed_rban_records = [
             Parsed_rBAN_Record.from_dict(record)
             for record in yaml.safe_load(f)
         ]

     nrp_variants = rban_records_to_nrp_variants(parsed_rban_records)
     return {
         variant.nrp_variant_id.nrp_id: variant
         for variant in nrp_variants
     }


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
    assert (nerpa_dir / 'nerpa.py').exists(), f"nerpa.py not found in {nerpa_dir}. Please ensure the script is run from the correct directory."

    args = load_command_line_args(nerpa_dir)

    logger = configure_logging(args.output_dir)

    # 1. Load configs and helpers
    nerpa_cfg = load_config()
    monomer_names_helper = load_monomer_names_helper(nerpa_cfg.monomers_config,
                                                     nerpa_dir)
    hmm_helper = load_hmm_helper(nerpa_cfg, monomer_names_helper)
    pnrpdb_info = pl.read_csv(args.pnrpdb2_info, sep='\t')

    # 2. Load approved matches and corresponding BGC variants
    logger.info('Loading approved matches')
    approved_matches = [
        TestMatch.from_yaml_dict(test_match_dict)
        for test_match_dict in yaml.safe_load(args.approved_matches.read_text())
    ]
    approved_matches = remove_deprecated_nrps(approved_matches, pnrpdb_info)

    logger.info('Loading BGC variants')
    test_matches_genome_ids = {
        m.bgc_id
        for m in approved_matches
    }
    bgc_variants = get_bgc_variants(
        args.nerpa_results_on_mibig_no_calibration,
        test_matches_genome_ids,
    )

    bgc_id_to_bgc_variants = defaultdict(list)
    for bgc_variant in bgc_variants:
        bgc_id_to_bgc_variants[bgc_variant.bgc_variant_id.bgc_id].append(bgc_variant)
    matched_bgc_variants: Dict[NRP_ID, BGC_Variant] = load_bgc_variants_for_matches(
        bgc_id_to_bgc_variants,
        approved_matches,
        logger
    )
    check_bgcs_with_many_variants(list(matched_bgc_variants.values()))

    nrp_id_to_nrp_variant = get_nrp_variants(args.mibig_norine_preprocessed)


    # note: alignments are taken from "old" approved matches while bgc predictions are taken from current matches
    matches_with_bgcs_nrps_for_training = [
        MatchWithBGCNRP(
            m,
            matched_bgc_variants[m.nrp_id],
            nrp_id_to_nrp_variant[m.nrp_id]
        )
        for m in approved_matches
        if m.nrp_id in matched_bgc_variants
    ]

    logger.info(f'Number of matches with BGC and NRP variants for training: {len(matches_with_bgcs_nrps_for_training)}')

    data_for_training = extract_data_for_training(
        matches_with_bgcs_nrps_for_training,
        hmm_helper,
        dir_for_hmm_figures=args.output_dir if args.draw_hmms else None
    )

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
