from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from src.matching.matcher_viterbi_algorithm import DetailedHMM
from typing import List, Tuple, Optional
from src.matching.matching_types_alignment import Alignment
from src.matching.matcher_viterbi_types import DetailedHMMStateType, DetailedHMMEdgeType
from src.generic.graphs import shortest_path_through
import networkx as nx
from src.data_types import NRP_Monomer
from src.matching.matcher_viterbi_types import DetailedHMMStateType, DetailedHMMEdgeType
from src.matching.matching_types_alignment import Alignment
from src.matching.matching_types_alignment_step import (
    AlignmentStep,
    AlignmentStep_BGC_Module_Info,
    AlignmentStep_NRP_Monomer_Info
)
from itertools import pairwise
from copy import deepcopy


def check_compatability(alignment: Alignment, hmm: DetailedHMM) -> bool:
    alignment_modules = [(step.bgc_module_info.gene_id, step.bgc_module_info.a_domain_idx)
                         for step in alignment
                         if step.bgc_module_info is not None]
    hmm_modules = [(module.gene_id, module.a_domain_idx)
                   for module in hmm.bgc_variant.modules]
    return alignment_modules == hmm_modules


def alignment_to_hmm_path(hmm: DetailedHMM, alignment: Alignment) -> List[Tuple[int, Optional[NRP_Monomer]]]:  # [(state_idx, emitted_monomer)]
    '''
    Given an alignment, return the path through the HMM that corresponds to it.
    The path is a sequence of state indices in the HMM.
    As I use it FOR training, edge weights are yet arbitrary and I don't use them.
    Instead, I find the shortest path by the number of edges,
    so that the path skippes the whole fragments or genes instead of individual modules.
    Also, I assign weight 0 to auxiliary edges so that only edges corresponding to the alignment are relevant.
    '''
    assert check_compatability(alignment, hmm), 'Alignment and HMM are incompatible'

    # I can't use step.step_type here, because of backwards compatibility issues
    match_steps_alignment_idxs = [i for i, step in enumerate(alignment)
                                  if step.bgc_module_info is not None and step.nrp_monomer_info is not None]
    match_steps = [alignment[i] for i in match_steps_alignment_idxs]

    # it could be that the same module is matched multiple times in case of iterations, I need all occurrences
    matched_module_idxs = [step.bgc_module_info.a_domain_idx for step in match_steps]

    module_match_state_idxs = [next(edge.to
                                    for edge in hmm.adj_list[hmm._module_idx_to_state_idx[module_idx]]
                                    if hmm.states[edge.to].state_type == DetailedHMMStateType.MATCH)
                               for module_idx in range(len(hmm.bgc_variant.modules))]

    # subpaths between adjacent matches
    paths_endpoints = [(hmm.start_state_idx,  # path before the first match
                        hmm._module_idx_to_state_idx[matched_module_idxs[0]])]
    for module_idx, next_module_idx in pairwise(matched_module_idxs):  # paths between matches
        paths_endpoints.append((module_match_state_idxs[module_idx],
                                hmm._module_idx_to_state_idx[next_module_idx]))
    paths_endpoints.append((module_match_state_idxs[-1],
                            hmm.final_state_idx))  # path after the last match

    # monomers inserted between matches
    monomers_subseqs = [alignment[step_idx].nrp_monomer_info  # monomers inserted before the first match
                        for step_idx in range(match_steps_alignment_idxs[0])]
    for step_idx, next_step_idx in pairwise(match_steps_alignment_idxs):
        monomers_subseqs.extend([alignment[i].nrp_monomer_info
                                 for i in range(step_idx, next_step_idx)
                                 if alignment[i].nrp_monomer_info is not None])
    monomers_subseqs.extend([alignment[step_idx].nrp_monomer_info  # monomers inserted after the last match
                             for step_idx in range(match_steps_alignment_idxs[-1], len(alignment))
                             if alignment[step_idx].nrp_monomer_info is not None])

    # edges which correspond to aligning steps with bgc_info != None or nrp_info != None
    main_edge_types = (
        DetailedHMMEdgeType.SKIP_FRAGMENT_AT_START,
        DetailedHMMEdgeType.INSERT_AT_START,
        DetailedHMMEdgeType.MATCH,
        DetailedHMMEdgeType.INSERT,
        DetailedHMMEdgeType.SKIP_MODULE,
        DetailedHMMEdgeType.SKIP_GENE,
        DetailedHMMEdgeType.SKIP_FRAGMENT,
        DetailedHMMEdgeType.SKIP_FRAGMENT_AT_END
    )


    path_with_emissions = []
    for (start, finish), emitted_monomers in zip(paths_endpoints, monomers_subseqs):
        # build path between consequent matches
        sub_hmm = deepcopy(hmm)
        for state_idx, state in enumerate(sub_hmm.states):
            if state.state_type == DetailedHMMStateType.MATCH and state_idx not in (start, finish):
                sub_hmm.adj_list[state_idx] = []  # by design there shouldn't be any matches in between so I make match states "deadends"
            # essentially, there's only one path between start and finish.
            # However, I need to make sure that the path skips the whole fragments or genes instead of individual modules
            # and chooses SKIP_FRAGMENT_AT_START/END instead of just SKIP_FRAGMENT when possible
            # that's why I assign these edge weights. Yes, dirty hacks :sweat_smile:
            for edge in sub_hmm.adj_list[state_idx]:
                match edge.edge_type:
                    case DetailedHMMEdgeType.SKIP_FRAGMENT_AT_START | DetailedHMMEdgeType.SKIP_FRAGMENT_AT_END:
                        edge.log_prob = 0
                    case DetailedHMMEdgeType.SKIP_FRAGMENT:
                        edge.log_prob = -1
                    case DetailedHMMEdgeType.SKIP_GENE:
                        edge.log_prob = -2
                    case DetailedHMMEdgeType.SKIP_MODULE:
                        edge.log_prob = -3
                    case _:
                        edge.log_prob = 0

        path_with_emissions.extend(sub_hmm.get_opt_path_with_emissions(start, finish, emitted_monomers))

    return path_with_emissions

