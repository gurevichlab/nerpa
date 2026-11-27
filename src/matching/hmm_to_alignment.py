from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from src.hmm.detailed_hmm import DetailedHMM
from typing import List
from src.hmm.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    DetailedHMMStateType
)
from src.hmm.hmm_scoring_types import ScoringSchemeName
from itertools import pairwise
from src.matching.alignment_type import Alignment
from src.matching.alignment_step_type import (
    AlignmentStep,
    AlignmentStep_BGC_Module_Info,
)
from src.rban_parsing.rban_monomer import rBAN_Monomer
from pathlib import Path

ET = DetailedHMMEdgeType  # for convenience
ST = DetailedHMMStateType

def hmm_path_to_alignment(
        hmm: DetailedHMM,
        path: List[int],
        nrp_monomers: List[rBAN_Monomer]
) -> Alignment:
    # for debugging
    # hmm.draw(Path('debug_hmm.svg'),
    #          highlight_path=path,
    #          emission_names=[
    #              mon.rban_name
    #              for mon in nrp_monomers
    #          ])

    scored_hmm = hmm.to_hmm(ScoringSchemeName.NERPA)
    monomers_iter = iter(nrp_monomers)
    alignment = []
    for edge_from, edge_to in pairwise(path):
        edge = hmm.transitions[edge_from][edge_to]

        # special case because multiple modules are affected
        if edge.edge_type in (
            ET.SKIP_UNTIL_NEXT_TENTATIVE_ASSEMBLY_LINE_START,
            ET.SKIP_MODULES_AT_END,
        ):
            fst_skipped_module_idx = (
                hmm.states[edge_from].related_module_idx
            )
            assert fst_skipped_module_idx is not None

            if (
                hmm.states[edge_to].state_type
                != ST.CHOOSE_IF_ITERATE_BGC
            ):
                next_module_idx = hmm.states[edge_to].related_module_idx
                assert next_module_idx is not None
                lst_skipped_module_idx = next_module_idx - 1
            else:
                lst_skipped_module_idx = (
                    len(hmm.bgc_variant.modules) - 1
                )

            alignment.extend(
                AlignmentStep(
                    bgc_module=AlignmentStep_BGC_Module_Info.from_bgc_module(hmm.bgc_variant.modules[module_idx]),
                    nrp_monomer=None,
                    score=edge.weight,
                    match_detailed_score=None,
                    step_type=edge.edge_type
                )
                for module_idx in range(fst_skipped_module_idx, lst_skipped_module_idx + 1)
            )
            continue


        if scored_hmm.emissions[edge_to]:
            emitted_monomer = next(monomers_iter)
            monomer_code = hmm.hmm_helper.monomer_names_helper.mon_to_int[
                emitted_monomer.to_base_mon()
            ]
            emission_score = scored_hmm.emissions[edge_to][monomer_code]
        else:
            emitted_monomer = None
            emission_score = 0

        # Note: a BGC module may be present even when this is not a match step.
        alignment_bgc_module_info = None
        match_detailed_score = None

        # Matching states: for these match_detailed_score is defined
        if hmm.states[edge_to].state_type in (
            ST.MATCH,
            ST.MATCH_POSSIBLE_ASSEMBLY_FINISH,
            ST.MATCH_LAST_MODULE
        ):
            module_idx = hmm.states[edge_to].related_module_idx
            assert module_idx is not None
            assert emitted_monomer is not None

            bgc_module = hmm.bgc_variant.modules[module_idx]
            alignment_bgc_module_info = (
                AlignmentStep_BGC_Module_Info.from_bgc_module(
                    bgc_module
                )
            )
            match_detailed_score = (
                hmm.hmm_helper.match_log_odds_components(
                    bgc_module,
                    emitted_monomer.to_base_mon(),
                    hmm.bgc_variant.has_pks_domains(),
                )
            )

        if edge.edge_type == ET.SKIP_MODULE:
            # SKIP_MODULE belongs to the module at edge_from.
            module_idx = hmm.states[edge_from].related_module_idx
            assert module_idx is not None
            bgc_module = hmm.bgc_variant.modules[module_idx]
            alignment_bgc_module_info = (
                AlignmentStep_BGC_Module_Info.from_bgc_module(
                    bgc_module
                )
            )

        alignment.append(
            AlignmentStep(
                bgc_module=alignment_bgc_module_info,
                nrp_monomer=emitted_monomer,
                score=edge.weight + emission_score,
                match_detailed_score=match_detailed_score,
                step_type=edge.edge_type
            )
        )

    return alignment
