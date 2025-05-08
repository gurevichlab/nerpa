from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from src.matching.detailed_hmm import DetailedHMM
from typing import List, Tuple, Optional
from src.matching.alignment_type import Alignment, show_alignment
from src.matching.hmm_auxiliary_types import DetailedHMMStateType, DetailedHMMEdgeType
from src.generic.graphs import shortest_path_through
import networkx as nx
from src.data_types import GeneId
from src.monomer_names_helper import NRP_Monomer
from src.matching.hmm_auxiliary_types import DetailedHMMStateType, DetailedHMMEdgeType
from src.matching.alignment_type import Alignment
from itertools import pairwise
from copy import deepcopy
from pathlib import Path


def alignment_to_hmm_path(hmm: DetailedHMM, alignment: Alignment) -> List[Tuple[int, Optional[NRP_Monomer]]]:  # [(state_idx, emitted_monomer)]

    '''
    Given an alignment, return the path through the HMM that corresponds to it.
    The path is a sequence of state indices in the HMM.
    As I use it FOR training, edge weights are yet arbitrary and I don't use them.
    Instead, I find the shortest path by the number of edges,
    so that the path skippes the whole fragments or genes instead of individual modules.
    Also, I assign weight 0 to auxiliary edges so that only edges corresponding to the alignment are relevant.
    '''
    # I can't use step.step_type here, because of backwards compatibility issues
    #print(f'Reconstructing hmm path for {hmm.bgc_variant.bgc_variant_id.bgc_id.genome_id}')
    #print(f'Alignment: \n{show_alignment(alignment)}')
    #hmm.draw(Path(f'hmm_{hmm.bgc_variant.bgc_variant_id.bgc_id.genome_id}.png'))
    match_steps_alignment_idxs = [i for i, step in enumerate(alignment)
                                  if step.bgc_module is not None and step.nrp_monomer is not None]
    match_steps = [alignment[i] for i in match_steps_alignment_idxs]

    # it could be that the same module is matched multiple times in case of iterations, I need all occurrences
    module_idx_in_bgc = {(module.gene_id, module.a_domain_idx): i
                         for i, module in enumerate(hmm.bgc_variant.modules)}
    matched_module_idxs = [module_idx_in_bgc[(step.bgc_module.gene_id, step.bgc_module.a_domain_idx)]
                           for step in match_steps]

    module_match_state_idxs = [next(edge_to
                                    for edge_to in hmm.transitions[hmm._module_idx_to_state_idx[module_idx]]
                                    if hmm.states[edge_to].state_type == DetailedHMMStateType.MATCH)
                               for module_idx in range(len(hmm.bgc_variant.modules))]

    # subpaths between adjacent matches
    paths_endpoints = [(hmm.start_state_idx,  # path before the first match
                        hmm._module_idx_to_state_idx[matched_module_idxs[0]])]
    for module_idx, next_module_idx in pairwise(matched_module_idxs):  # paths between matches
        paths_endpoints.append((module_match_state_idxs[module_idx],
                                hmm._module_idx_to_state_idx[next_module_idx]))
    paths_endpoints.append((module_match_state_idxs[matched_module_idxs[-1]],
                            hmm.final_state_idx))  # path after the last match

    # monomers inserted between matches
    monomers_subseqs = [[alignment[step_idx].nrp_monomer  # monomers inserted before the first match
                        for step_idx in range(match_steps_alignment_idxs[0])
                         if alignment[step_idx].nrp_monomer is not None]]
    for step_idx, next_step_idx in pairwise(match_steps_alignment_idxs):
        monomers_subseqs.append([alignment[i].nrp_monomer
                                 for i in range(step_idx, next_step_idx)
                                 if alignment[i].nrp_monomer is not None])
    monomers_subseqs.append([alignment[step_idx].nrp_monomer  # monomers inserted after the last match
                             for step_idx in range(match_steps_alignment_idxs[-1], len(alignment))
                             if alignment[step_idx].nrp_monomer is not None])

    # edges which correspond to aligning steps with bgc_info != None or nrp_info != None

    path_with_emissions = []
    for (start, finish), emitted_monomers in zip(paths_endpoints, monomers_subseqs):
        #print(start, finish, emitted_monomers)
        # build path between consequent matches
        sub_hmm = deepcopy(hmm)
        for state in sub_hmm.states:
            for mon in state.emissions:
                state.emissions[mon] = 0

        for state_idx, state in enumerate(sub_hmm.states):
            if state.state_type == DetailedHMMStateType.MATCH and state_idx != start:
                sub_hmm.transitions[state_idx] = {}  # by design there shouldn't be any matches in between so I make match states "deadends"
            # essentially, there's only one path between start and finish.
            # However, I need to make sure that the path skips the whole fragments or genes instead of individual modules
            # and chooses SKIP_FRAGMENT_AT_START/END instead of just SKIP_FRAGMENT when possible
            # that's why I assign these edge weights. Yes, dirty hacks :sweat_smile:
            for edge_to, edge_info in sub_hmm.transitions[state_idx].items():
                match edge_info.edge_type:
                    case DetailedHMMEdgeType.SKIP_GENE:
                        log_prob = -2
                    case DetailedHMMEdgeType.SKIP_MODULE:
                        log_prob = -3
                    case _:
                        log_prob = 0
                sub_hmm.transitions[state_idx][edge_to] = sub_hmm.transitions[state_idx][edge_to]._replace(weight=log_prob)
        sub_hmm._hmm = None
        try:
            path_with_emissions.extend(sub_hmm.get_opt_path_with_emissions(start, finish, emitted_monomers))
        except:
            print(f'start: {start}, finish: {finish}, emitted_monomers: {emitted_monomers}')
            raise

    return path_with_emissions

