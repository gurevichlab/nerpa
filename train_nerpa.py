from __future__ import annotations
import argparse
from argparse import Namespace as CommandlineArgs
from collections import defaultdict
from logging import Logger
from pathlib import Path
import traceback
from typing import Dict, List, NamedTuple, Set, Optional

import polars as pl
import yaml
import json

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
from src.hmm.hmm_scoring_helper import HMMHelper
from src.matching.matcher import NRP_ID
from src.monomer_names_helper.monomer_names_helper import MonomerNamesHelper
from src.rban_parsing.nrp_variant_types import NRP_Variant
from src.rban_parsing.rban_parser import Parsed_rBAN_Record
from src.rban_parsing.retrieve_nrp_variants import rban_records_to_nrp_variants
from src.testing.testing_types import (
    CuratedAlignment,
    CuratedAlignmentId,
)
from src.training.hmm_parameters.bgc_alignment_compatibility import (
    bgc_variant_match_compatible,
)
from src.training.hmm_parameters.extract_data_for_training import (
    extract_data_for_training,
)
from src.training.hmm_parameters.transition_probabilities.estimate_transition_probs import estimate_hmm_transition_probs
from src.training.hmm_parameters.emission_probabilities.hmm_infer_emission_params import infer_emission_params
from src.training.hmm_parameters.norine_monomers_info import NorineMonomersInfo
from src.training.hmm_parameters.training_types import AlignmentWithBGCNRP
from src.training.write_results import write_params
from src.training.logging_config import configure_logging
from src.training.hmm_parameters.transition_probabilities.estimate_event_prob_given_context import EstimateEventProbGivenContextConfig
from src.training.hmm_parameters.emission_probabilities.hmm_infer_emission_params import StepFunctionEstimationConfig
from src.training.hmm_parameters.norine_monomers_info import NorineMonomersInfo
from src.training.hmm_parameters.short_aliases import *
from src.training.hmm_parameters.emission_probabilities.hmm_infer_emission_params import EmissionParams
from src.training.hmm_parameters.training_types import (
    HMM_Chosen_Steps_df,
    Estimated_Transition_Probs,
)
from src.general_type_aliases import LogProb, Prob
from src.hmm.hmm_scoring_config import HMMScoringConfig


MIBIG_ID = str

def load_bgc_variants_by_mibig_id(
        nerpa_results: Path,
        mibig_ids_to_include: Optional[Set[str]]=None,
) -> Dict[MIBIG_ID, List[BGC_Variant]]:
    if mibig_ids_to_include is None:
        mibig_ids_to_include = set()
        
    bgc_variants_yaml = nerpa_results / 'preprocessed_input' / 'BGC_variants.yaml'
    bgc_variants = defaultdict(list)
    yaml_variants = [
        BGC_Variant.from_dict(bgc_variant_dict)
        for bgc_variant_dict in yaml.safe_load(bgc_variants_yaml.read_text())
    ]
    for bgc_variant in yaml_variants:
        genome_id = bgc_variant.bgc_variant_id.bgc_id._genome_id()
        if genome_id in mibig_ids_to_include:
            bgc_variants[genome_id].append(bgc_variant)
    return bgc_variants


def get_alignment_id_to_bgc_variant(
        mibig_id_to_bgc_variants: Dict[MIBIG_ID, List[BGC_Variant]],
        curated_alignments: List[CuratedAlignment],
        log: Logger
) -> Dict[CuratedAlignmentId, BGC_Variant]:
    alignment_id_to_bgc_variant = {}
    for m in curated_alignments:
        try:
            bgc_variant = next(
                bgc_variant
                for bgc_variant in mibig_id_to_bgc_variants[m.mibig_id]
                if bgc_variant_match_compatible(bgc_variant, m)
            )
        except StopIteration:
            log.warning(f'Cannot find appropriate BGC variant for {m.nrp_id}')
            log.debug(f'Match that failed to find variant: {m}')
            continue

        alignment_id_to_bgc_variant[m.get_id()] = bgc_variant

    return alignment_id_to_bgc_variant


def get_nrp_id_to_nrp_variant(parsed_rban_records_yaml: Path) -> Dict[NRP_ID, NRP_Variant]:
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


def load_hmm_helper(
        cfg: Config,
        monomer_names_helper: MonomerNamesHelper,
        norine_monomers_info: NorineMonomersInfo
) -> HMMHelper:
    specificity_prediction_helper = SpecificityPredictionHelper(
        cfg.specificity_prediction_config,
        monomer_names_helper
    )

    hmm_scoring_config = HMMScoringConfig(
        nerpa_dir=cfg.nerpa_dir,
        path_to_config=cfg.hmm_scoring_config,
        specificity_prediction_helper=specificity_prediction_helper,
        norine_monomers_info=norine_monomers_info
    )

    return HMMHelper(hmm_scoring_config, monomer_names_helper)


def parse_cli_args(nerpa_dir: Path) -> CommandlineArgs:
    parser = argparse.ArgumentParser(description="Estimates Nerpa parameters "
                                                 "based on Nerpa results on approved matches")
    default_curated_alignments = (
        nerpa_dir
        / 'data'
        / 'for_training_and_testing'
        / 'curated_alignments.json'
    )
    parser.add_argument(
        "--curated-alignments",
        type=Path,
        default=default_curated_alignments,
        help=f"Path to the JSON with curated alignments file (default: {default_curated_alignments})",
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
    parser.add_argument(
        "--norine-monomers-info",
        type=Path,
        help=f"Path to the TSV file containing Norine monomers info",
    )
        
    parser.add_argument("--output-dir", type=Path,
                        default=nerpa_dir / 'training_results')
    parser.add_argument("--draw-hmms", action='store_true')
    return parser.parse_args()


def check_bgcs_with_many_variants(bgc_variants: List[BGC_Variant]):
    bgc_id_to_variants = defaultdict(set)
    for bgc_variant in bgc_variants:
        bgc_id_to_variants[bgc_variant.bgc_variant_id.bgc_id].add(bgc_variant.bgc_variant_id)

    for bgc_id, variant_ids in bgc_id_to_variants.items():
        if len(variant_ids) > 1:
            print(
                f'Warning: BGC {bgc_id} has {len(variant_ids)} '
                f'variants used for training:\n {variant_ids}\n'
            )


def main():
    nerpa_dir = Path(__file__).parent
    assert (nerpa_dir / 'nerpa.py').exists(), f"nerpa.py not found in {nerpa_dir}. Please ensure the script is run from the correct directory."

    args = parse_cli_args(nerpa_dir)

    logger = configure_logging(args.output_dir / 'train_nerpa.log')

    # 1. Load configs and helpers
    nerpa_cfg = load_config()
    monomer_names_helper = load_monomer_names_helper(nerpa_cfg.monomers_config, nerpa_dir)
    norine_monomers_info: NorineMonomersInfo = (
       NorineMonomersInfo.from_tsv(args.norine_monomers_info)
    )
    hmm_helper = load_hmm_helper(
        nerpa_cfg,
        monomer_names_helper,
        norine_monomers_info
    )
    pnrpdb_info = pl.read_csv(args.pnrpdb2_info, separator='\t')

    # 2. Load curated alignments and corresponding BGC variants
    logger.info('Loading curated alignments')
    curated_alignments_data = json.loads(args.curated_alignments.read_text(encoding="utf-8"))
    curated_alignments: List[CuratedAlignment] = [
        CuratedAlignment.from_yaml_dict(ca)
        for ca in curated_alignments_data["curated_alignments"]
    ]
    # curated_alignments = remove_deprecated_nrps(curated_alignments, pnrpdb_info)

    logger.info('Loading BGC variants')
    test_matches_genome_ids = {
        m.mibig_id
        for m in curated_alignments
    }
    mibig_id_to_bgc_variants = load_bgc_variants_by_mibig_id(
        args.nerpa_results_on_mibig_no_calibration,
        test_matches_genome_ids,
    )

    alignment_id_to_bgc_variant: Dict[CuratedAlignmentId, BGC_Variant] = get_alignment_id_to_bgc_variant(
        mibig_id_to_bgc_variants,
        curated_alignments,
        logger
    )
    check_bgcs_with_many_variants(list(alignment_id_to_bgc_variant.values()))

    nrp_id_to_nrp_variant = get_nrp_id_to_nrp_variant(args.mibig_norine_preprocessed)

    # note: alignments are taken from "old" approved matches while bgc predictions are taken from current matches
    alignments_with_bgcs_nrps_for_training = [
        AlignmentWithBGCNRP(
            al,
            alignment_id_to_bgc_variant[al.get_id()],
            nrp_id_to_nrp_variant[al.nrp_id]
        )
        for al in curated_alignments
        if al.get_id() in alignment_id_to_bgc_variant
        and al.nrp_id in nrp_id_to_nrp_variant
    ]

    logger.info(f'Number of matches with BGC and NRP variants for training: {len(alignments_with_bgcs_nrps_for_training)}')

    hmm_chosen_steps: HMM_Chosen_Steps_df = extract_data_for_training(
        alignments_with_bgcs_nrps_for_training,
        hmm_helper,
        logger,
        dir_for_hmm_figures=(
            args.output_dir
            if args.draw_hmms
            else None
        )
    )

    # Dump extracted data for debug
    hmm_chosen_steps.write_tsv(args.output_dir / 'hmm_chosen_steps.tsv')
    hmm_chosen_steps.write_transition_cnts_tsv(args.output_dir / 'hmm_transition_counts.tsv')


    hmm_transition_probs_cfg = EstimateEventProbGivenContextConfig()
    estimated_transition_probs: Estimated_Transition_Probs = (
        estimate_hmm_transition_probs(
            hmm_chosen_steps,
            cfg=hmm_transition_probs_cfg,
            logger=logger,
        )
    )

    emission_params: EmissionParams = infer_emission_params(
        hmm_chosen_steps_df=hmm_chosen_steps,
        norine_monomers_info=norine_monomers_info,
        step_function_cfg=StepFunctionEstimationConfig(),
        logger=logger,
        output_dir=args.output_dir,   # for plotting step function
    )
    write_params(
        estimated_transition_probs=estimated_transition_probs,
        emission_params=emission_params,
        data_for_training=hmm_chosen_steps,
        old_configs_dir=nerpa_dir / 'configs',
        output_dir=args.output_dir
    )


if __name__ == "__main__":
    try:
        main()
    except Exception:
        with Path("train_nerpa_err.log").open("w") as error_file:
            traceback.print_exc(file=error_file)
        raise

