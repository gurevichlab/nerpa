from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from src.matching.matcher_viterbi_algorithm import DetailedHMM
from typing import List
from src.matching.matcher_viterbi_types import (
    DetailedHMMEdgeType,
    DetailedHMMStateType,
    DetailedHMMState,
    DetailedHMMEdge,
)
from src.matching.matching_types_alignment_step import AlignmentStep
from itertools import pairwise
from src.data_types import NRP_Monomer
from src.matching.matching_types_alignment import Alignment
from src.matching.matching_types_alignment_step import (
    AlignmentStep,
    AlignmentStep_BGC_Module_Info,
    AlignmentStep_NRP_Monomer_Info
)
from src.matching.scoring_helper import ScoringHelper


def hmm_path_to_alignment(hmm: DetailedHMM,
                          path: List[int],
                          nrp_monomers: List[NRP_Monomer],
                          scoring_helper: ScoringHelper) -> Alignment:
    # The idea is that first the transition is made and then the emission -- in the new state
    # TODO: maybe refactor to do first emission and then transition to make it more similar to the standard Viterbi algorithm
    alignment = []
    mon_idx = 0
    for edge_from, edge_to in pairwise(path):
        edge = next(edge for edge in hmm.adj_list[edge_from] if edge.to == edge_to)  # maybe use dicts for speedup?

        # these edge types are handled separately because they involve multiple modules
        # and hence yield multiple alignment steps
        if edge.edge_type in (DetailedHMMEdgeType.SKIP_GENE,
                              DetailedHMMEdgeType.SKIP_FRAGMENT,
                              DetailedHMMEdgeType.SKIP_FRAGMENT_AT_START,
                              DetailedHMMEdgeType.SKIP_FRAGMENT_AT_END):
            fst_module_idx = hmm.state_idx_to_module_idx[edge_from]
            last_module_idx = hmm.state_idx_to_module_idx[edge_to]
            alignment.extend(AlignmentStep(bgc_module_info=AlignmentStep_BGC_Module_Info(hmm.bgc_variant.modules[module_idx]),
                                           nrp_monomer_info=None,
                                           score=edge.log_prob,
                                           match_detailed_score=None,
                                           step_type=edge.edge_type)
                             for module_idx in range(fst_module_idx, last_module_idx + 1))
            continue

        bgc_module = None
        nrp_monomer = None
        match_detailed_score = None
        transition_score = edge.log_prob
        emission_score = 0

        match edge.edge_type:
            # with emissions
            case DetailedHMMEdgeType.MATCH:
                bgc_module = hmm.bgc_variant.modules[hmm.state_idx_to_module_idx[edge_from]]
                nrp_monomer = nrp_monomers[mon_idx]
                match_detailed_score = scoring_helper.match_detailed_score(bgc_module, nrp_monomers[mon_idx])
                emission_score = sum(match_detailed_score)
                mon_idx += 1
            case DetailedHMMEdgeType.START_INSERTING | DetailedHMMEdgeType.INSERT:
                bgc_module = hmm.bgc_variant.modules[hmm.state_idx_to_module_idx[edge_from]]
                nrp_monomer = nrp_monomers[mon_idx]
                # in the current setup, emission_score is always 0 for insertions, as if the emission is always the same
                mon_idx += 1
            case DetailedHMMEdgeType.START_INSERTING_AT_START | DetailedHMMEdgeType.INSERT_AT_START:
                next_module_state_idx = next(next_edge.to
                                             for next_edge in hmm.adj_list[edge_to]
                                             if hmm.states[next_edge.to].state_type == DetailedHMMStateType.MODULE_START)
                bgc_module = hmm.bgc_variant.modules[next_module_state_idx]
                # in the current setup, emission_score is always 0 for insertions, as if the emission is always the same
                nrp_monomer = nrp_monomers[mon_idx]
                mon_idx += 1
            # skips
            case DetailedHMMEdgeType.SKIP_MODULE | DetailedHMMEdgeType.SKIP_GENE | DetailedHMMEdgeType.SKIP_FRAGMENT:
                bgc_module = hmm.bgc_variant.modules[hmm.state_idx_to_module_idx[edge_from]]

        score = transition_score + emission_score
        bgc_module_info = AlignmentStep_BGC_Module_Info.from_bgc_module(bgc_module) \
            if bgc_module is not None else None
        nrp_monomer_info = AlignmentStep_NRP_Monomer_Info.from_nrp_monomer(nrp_monomer) \
            if nrp_monomer is not None else None
        alignment.append(AlignmentStep(bgc_module_info=bgc_module_info,
                                       nrp_monomer_info=nrp_monomer_info,
                                       score=score,
                                       match_detailed_score=match_detailed_score,
                                       step_type=edge.edge_type))

    return alignment