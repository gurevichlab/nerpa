from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from src.matching.matcher_viterbi_detailed_hmm import DetailedHMM
from typing import List, Tuple, Optional
from src.matching.matching_types_alignment import Alignment
from src.matching.matcher_viterbi_types import DetailedHMMStateType, DetailedHMMEdgeType
from src.generic.graphs import shortest_path_through
import networkx as nx
from src.data_types import GeneId
from src.monomer_names_helper import NRP_Monomer, MonomerNamesHelper, monomer_names_helper
from src.matching.matcher_viterbi_types import DetailedHMMStateType, DetailedHMMEdgeType
from src.matching.matching_types_alignment import Alignment
from itertools import pairwise
from copy import deepcopy
from pathlib import Path


def fix_indexing_in_alignment(alignment: Alignment):
    if min(step.bgc_module.a_domain_idx
           for step in alignment
           if step.bgc_module is not None) == 1:
        for step in alignment:
            if step.bgc_module is not None:
                step.bgc_module = step.bgc_module._replace(a_domain_idx=step.bgc_module.a_domain_idx - 1)

def fix_gene_names_in_alignment(alignment: Alignment):
    if all(step.bgc_module.gene_id.startswith('ctg')
           for step in alignment
           if step.bgc_module is not None) == 1:
        for step in alignment:
            if step.bgc_module is not None:
                step.bgc_module = step.bgc_module._replace(gene_id=GeneId(step.bgc_module.gene_id.split('_')[1]))

def fix_residue_names_in_alignment(alignment: Alignment):
    for step in alignment:
        if step.nrp_monomer is not None:
            parsed_monomer = monomer_names_helper.parsed_name(step.nrp_monomer.rban_name, name_format='norine')
            step.nrp_monomer = step.nrp_monomer._replace(residue=parsed_monomer.residue)


def check_compatability(alignment: Alignment, hmm: DetailedHMM) -> bool:
    alignment_modules = [(step.bgc_module.gene_id, step.bgc_module.a_domain_idx)
                         for step in alignment
                         if step.bgc_module is not None]
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
    # I can't use step.step_type here, because of backwards compatibility issues
    print(hmm.bgc_variant.genome_id)
    match_steps_alignment_idxs = [i for i, step in enumerate(alignment)
                                  if step.bgc_module is not None and step.nrp_monomer is not None]
    match_steps = [alignment[i] for i in match_steps_alignment_idxs]

    # it could be that the same module is matched multiple times in case of iterations, I need all occurrences
    module_idx_in_bgc = {(module.gene_id, module.a_domain_idx): i
                         for i, module in enumerate(hmm.bgc_variant.modules)}
    matched_module_idxs = [module_idx_in_bgc[(step.bgc_module.gene_id, step.bgc_module.a_domain_idx)]
                           for step in match_steps]

    module_match_state_idxs = [next(edge_to
                                    for edge_to in hmm.adj_list[hmm._module_idx_to_state_idx[module_idx]]
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
    cnt = 0
    for (start, finish), emitted_monomers in zip(paths_endpoints, monomers_subseqs):
        # build path between consequent matches
        cnt += 1
        sub_hmm = deepcopy(hmm)
        for state_idx, state in enumerate(sub_hmm.states):
            if state.state_type == DetailedHMMStateType.MATCH and state_idx != start:
                sub_hmm.adj_list[state_idx] = {}  # by design there shouldn't be any matches in between so I make match states "deadends"
            # essentially, there's only one path between start and finish.
            # However, I need to make sure that the path skips the whole fragments or genes instead of individual modules
            # and chooses SKIP_FRAGMENT_AT_START/END instead of just SKIP_FRAGMENT when possible
            # that's why I assign these edge weights. Yes, dirty hacks :sweat_smile:
            for edge_to, edge_info in sub_hmm.adj_list[state_idx].items():
                match edge_info.edge_type:
                    case DetailedHMMEdgeType.SKIP_FRAGMENT_AT_START | DetailedHMMEdgeType.SKIP_FRAGMENT_AT_END:
                        log_prob = 0
                    case DetailedHMMEdgeType.SKIP_FRAGMENT:
                        log_prob = -1
                    case DetailedHMMEdgeType.SKIP_GENE:
                        log_prob = -2
                    case DetailedHMMEdgeType.SKIP_MODULE:
                        log_prob = -3
                    case _:
                        log_prob = 0
                sub_hmm.adj_list[state_idx][edge_to] = sub_hmm.adj_list[state_idx][edge_to]._replace(log_prob=log_prob)
        # sub_hmm.draw(Path(f'sub_hmm_{cnt}.png'))
        path_with_emissions.extend(sub_hmm.get_opt_path_with_emissions(start, finish, emitted_monomers))
        path_with_emissions.append((finish, None))  # add the last state without emitted monomer

    return path_with_emissions

