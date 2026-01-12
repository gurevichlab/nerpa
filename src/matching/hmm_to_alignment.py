from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from src.hmm.detailed_hmm import DetailedHMM
from typing import List
from src.hmm.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    DetailedHMMStateType
)
from itertools import pairwise
from src.matching.alignment_type import Alignment
from src.matching.alignment_step_type import (
    AlignmentStep,
    AlignmentStep_BGC_Module_Info,
)
from src.rban_parsing.rban_monomer import rBAN_Monomer

ET = DetailedHMMEdgeType  # for convenience
ST = DetailedHMMStateType

def hmm_path_to_alignment(hmm: DetailedHMM,
                          path: List[int],
                          nrp_monomers: List[rBAN_Monomer]) -> Alignment:
    # for u in path:
    #     print(f'State {u}: {hmm.states[u].state_type}, related module idx: {hmm.states[u].related_module_idx}')
    monomers_iter = iter(nrp_monomers)
    alignment = []
    for edge_from, edge_to in pairwise(path):
        edge = hmm.transitions[edge_from][edge_to]

        # special case because multiple modules are affected
        if edge.edge_type in (
            ET.SKIP_UNTIL_NEXT_TENTATIVE_ASSEMBLY_LINE_START,
            ET.SKIP_MODULES_AT_END,
        ):
            fst_skipped_module_idx = hmm.states[edge_from].related_module_idx
            lst_skipped_module_idx = (
                    hmm.states[edge_to].related_module_idx - 1
                    if hmm.states[edge_to].state_type != ST.FINAL
                    else len(hmm.bgc_variant.modules) - 1
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


        if hmm.states[edge_to].emissions:
            emitted_monomer = next(monomers_iter)
            emission_score = hmm.states[edge_to].emissions[emitted_monomer.to_base_mon()]
        else:
            emitted_monomer = None
            emission_score = 0

        # default values
        alignment_bgc_module_info = None  # Note: alignment_bgc_module_info can be None even if bgc_module is not None
        match_detailed_score = 0

        # Matching states: for these match_detailed_score is defined
        if hmm.states[edge_to].state_type in (
            ST.MATCH,
            ST.MATCH_POSSIBLE_ASSEMBLY_FINISH,
            ST.MATCH_LAST_MODULE
        ):
            bgc_module = hmm.bgc_variant.modules[hmm.states[edge_to].related_module_idx]
            alignment_bgc_module_info = AlignmentStep_BGC_Module_Info.from_bgc_module(bgc_module)
            match_detailed_score = hmm.hmm_helper.normalized_match_detailed_score(bgc_module,
                                                                                  emitted_monomer.to_base_mon(),
                                                                                  hmm.bgc_variant.has_pks_domains())

        if edge.edge_type == ET.SKIP_MODULE:
            bgc_module = hmm.bgc_variant.modules[hmm.states[edge_from].related_module_idx]  # careful! edge_from used here
            alignment_bgc_module_info = AlignmentStep_BGC_Module_Info.from_bgc_module(bgc_module)

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