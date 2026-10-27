from pathlib import Path
from typing import (
    Dict,
    List,
    Tuple,
    Optional,
    NamedTuple,
)

import pandas as pd
import polars as pl
import logging

from src.antismash_parsing.antismash_parser_types import BGC_ID
from src.matching.alignment_type import light_alignment_to_str
from src.hmm.hmm_auxiliary_types import (
    DetailedHMMStateType,
    DetailedHMMEdgeType,
    StateIdx,
)
from src.hmm.hmm_constructor.hmm_constructor_state_edge_context_relations import (
    MATCHING_STATE_TYPES,
    RELEVANT_GENOMIC_CONTEXT,
    INSERT_STATE_TYPES,
)
from src.hmm.hmm_scoring_helper import HMMHelper
from src.hmm.detailed_hmm import DetailedHMM
from src.rban_parsing.nrp_variant_types import NRP_Variant
from src.rban_parsing.rban_monomer import rBAN_Monomer
from src.testing.simplified_alignment import simplified_alignment_to_light_alignments, simplified_alignment_to_str
from src.training.hmm_parameters.training_types import (
    AlignmentWithBGCNRP,
    HMM_Chosen_Steps_df,
    HMM_Chosen_Steps_Schema,
)
from src.monomer_names_helper import NRP_Monomer
from itertools import pairwise
from collections import defaultdict

PathWithEmissions = List[Tuple[StateIdx, Optional[NRP_Monomer]]]

class HMM_Match_Data(NamedTuple):
    detailed_hmm: DetailedHMM
    nrp_variant: NRP_Variant
    path_with_emissions: PathWithEmissions


def get_hmm_steps(
        detailed_hmm: DetailedHMM,
        nrp_variant: NRP_Variant,
        path_with_emissions: PathWithEmissions
) -> HMM_Chosen_Steps_df:
    """
    At each turn of the path a choice between outcoming edges is made.
    This function returns the information about outcoming edges at each turn of the path,
    and indications whether the edge was chosen or not
    """
    ET = DetailedHMMEdgeType
    ECS = HMM_Chosen_Steps_Schema
    bgc_modules = detailed_hmm.bgc_variant.modules

    idx_in_insert_run: int | None = None
    edge_choices = []
    for (u, emission), (v, _) in pairwise(path_with_emissions):
        state = detailed_hmm.states[u]
        module_idx = state.related_module_idx
        module = (
            bgc_modules[module_idx]
            if module_idx is not None
            else None
        )
        edge = detailed_hmm.transitions[u][v]
        idx_in_insert_run = (
            ((idx_in_insert_run or 0) + 1)
            if edge.edge_type in INSERT_STATE_TYPES
            else None
        )

        edge_choices.append({
            ECS.BGC_ID: detailed_hmm.bgc_variant.bgc_variant_id,
            ECS.NRP_ID: nrp_variant.nrp_variant_id,
            ECS.STATE_TYPE: state.state_type,
            ECS.EDGE_TYPE: edge.edge_type,
            ECS.IDX_IN_INSERT_RUN: idx_in_insert_run,
            ECS.EDGE_FROM: u,
            ECS.EDGE_TO: v,
            ECS.MODULE_IDX: module_idx,
            ECS.MODULE: module,
            ECS.EMISSION: emission,
        })


    return HMM_Chosen_Steps_df.from_rows(edge_choices)


def get_all_hmm_chosen_steps(
        hmms_with_paths_with_emissions: List[HMM_Match_Data]
) -> HMM_Chosen_Steps_df:
    hmm_chosen_steps = pl.DataFrame()

    for nrp_variant, detailed_hmm, path_with_emissions in hmms_with_paths_with_emissions:
        new_chosen_steps = get_hmm_steps(
            nrp_variant,
            detailed_hmm,
            path_with_emissions
        ).df
        hmm_chosen_steps = pl.concat([hmm_chosen_steps, new_chosen_steps], how="vertical")

    return HMM_Chosen_Steps_df(hmm_chosen_steps)

    
def get_hmms_with_paths_with_emissions(
        alignments_with_bgcs_nrps: List[AlignmentWithBGCNRP],
        hmm_helper: HMMHelper,
        logger: logging.Logger,
        dir_for_hmm_figures: Optional[Path] = None
) -> List[HMM_Match_Data]:
    hmms_with_paths_with_emissions = []
    for al, bgc_variant, nrp_variant in alignments_with_bgcs_nrps:
        try:
            detailed_hmm = DetailedHMM.from_bgc_variant(bgc_variant, hmm_helper)
        except Exception as e:
            logger.warning(f'WARNING: cannot construct DetailedHMM for match {al.nrp_id} ')
            raise

        all_alignments = list(simplified_alignment_to_light_alignments(
            al.true_alignment,
            detailed_hmm.bgc_variant,
            nrp_variant
        ))
        if not all_alignments:
            logger.warning(
                f'WARNING: cannot reconstruct alignment for match {al.nrp_id}\n'
                'Original alignment:\n'
                f'{simplified_alignment_to_str(al.true_alignment)}\n'
                'NRP variant:\n'
                f'{nrp_variant.to_str_compact()}'
            )
            continue

        if len(all_alignments) > 1:
            first_alignment_str = "\n".join(
                light_alignment_to_str(light_al) for light_al in all_alignments[0]
            )
            other_alignment_str = "\n".join(
                light_alignment_to_str(light_al) for light_al in all_alignments[1]
            )
            logger.warning(
                "\n".join(
                    [
                        f"WARNING: multiple alignments for match {al.nrp_id}",
                        "True alignment:",
                        simplified_alignment_to_str(al.true_alignment),
                        "First reconstructed alignment:",
                        first_alignment_str,
                        "Other reconstructed alignment:",
                        other_alignment_str,
                        f"WARNING: ambiguous alignment reconstruction for match {al.nrp_id} picking the first one",
                    ]
                )
            )

        alignments = all_alignments[0]
        for i, alignment in enumerate(alignments):
            try:
                path_with_emissions = detailed_hmm.alignment_to_path_with_emisions(
                    alignment,
                    alignment_type='light'
                )

                if dir_for_hmm_figures is not None:
                    detailed_hmm.draw(dir_for_hmm_figures / Path(f'{detailed_hmm.bgc_variant.bgc_variant_id.bgc_id.antiSMASH_file}_path_{i}.png'),
                                      highlight_path=[state for state, _ in path_with_emissions])  # for debugging
            except Exception as e:
                logger.warning(
                    f"WARNING: alignment {i} for match {al.get_id()} "
                    f"could not be converted to path with emissions: {e}"
                )
                continue

            hmms_with_paths_with_emissions.append(HMM_Match_Data(
                detailed_hmm,
                nrp_variant,
                path_with_emissions
            ))

    return hmms_with_paths_with_emissions


def extract_data_for_training(
        alignments_with_bgcs_nrps: List[AlignmentWithBGCNRP],
        hmm_helper: HMMHelper,
        logger: logging.Logger,
        dir_for_hmm_figures: Optional[Path] = None
) -> HMM_Chosen_Steps_df:
    hmms_nrps_paths: List[HMM_Match_Data] = (
        get_hmms_with_paths_with_emissions(
            alignments_with_bgcs_nrps,
            hmm_helper,
            logger,
            dir_for_hmm_figures
        )
    )
    hmms_chosen_steps_df: HMM_Chosen_Steps_df = (
        get_all_hmm_chosen_steps(hmms_nrps_paths)
    )
    return hmms_chosen_steps_df



