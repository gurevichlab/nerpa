from __future__ import annotations
from typing import List, Tuple, NamedTuple, Dict, Optional
from dataclasses import dataclass
from enum import Enum, auto

from src.data_types import NRP_Monomer, BGC_Module, BGC_Variant, GeneId
from src.matching.matcher_viterbi_types import DetailedHMMEdgeType, DetailedHMMStateType, DetailedHMMState, DetailedHMMEdge, HMM
from src.matching.auxilary import get_genes_intervals, get_fragments_intervals, get_emissions
from collections import defaultdict


def add_modules(states: List[DetailedHMMState],
                adj_list: Dict[int, List[DetailedHMMEdge]],
                bgc_variant: BGC_Variant,
                final_state_idx: int = -1) -> Dict[int, int]:  # module_idx -> module_start_state_idx

    # MODULE STATES: MODULE_START -> MATCH -> INSERT_MONOMER -> MODULE_START -> ...
    for i, module in enumerate(bgc_variant.modules):
        # add states
        states.append(DetailedHMMState(state_type=DetailedHMMStateType.MODULE_START,
                                       emissions={}))
        states.append(DetailedHMMState(state_type=DetailedHMMStateType.MATCH,
                                       emissions=get_emissions(module)))
        states.append(DetailedHMMState(state_type=DetailedHMMStateType.INSERT_MONOMER,
                                       emissions={}))  # TODO: fill in

        # add edges
        start_idx, match_idx, insert_idx = 3 * i + 1, 3 * i + 2, 3 * i + 3
        next_module_idx = 3 * i + 4 if i + 1 < len(bgc_variant.modules) else final_state_idx
        adj_list[start_idx] = [DetailedHMMEdge(edge_type=DetailedHMMEdgeType.MATCH,
                                               log_prob=0,  # TODO: fill in
                                               to=match_idx),
                               DetailedHMMEdge(edge_type=DetailedHMMEdgeType.SKIP_MODULE,
                                               log_prob=0,  # TODO: fill in
                                               to=next_module_idx)]
        adj_list[match_idx] = [DetailedHMMEdge(edge_type=DetailedHMMEdgeType.NO_INSERTIONS,
                                               log_prob=0,  # TODO: fill in
                                               to=next_module_idx),
                               DetailedHMMEdge(edge_type=DetailedHMMEdgeType.START_INSERTING,
                                               log_prob=0,  # TODO: fill in
                                               to=insert_idx)]
        adj_list[insert_idx] = [DetailedHMMEdge(edge_type=DetailedHMMEdgeType.CONTINUE_INSERTING,
                                                log_prob=0,  # TODO: fill in
                                                to=insert_idx),
                                DetailedHMMEdge(edge_type=DetailedHMMEdgeType.END_INSERTING,
                                                log_prob=0,  # TODO: fill in
                                                to=next_module_idx)]
    return {i: 3 * i + 1 for i in range(len(bgc_variant.modules))}


def add_skip_gene_and_fragments_middle(adj_list: Dict[int, List[DetailedHMMEdge]],
                                       bgc_variant: BGC_Variant,
                                       gene_intervals: Dict[GeneId, Tuple[int, int]],
                                       fragment_intervals: Dict[int, Tuple[int, int]],
                                       module_idx_to_start_state_idx: Dict[int, int]) -> None:
    # SKIP GENES AND FRAGMENTS IN THE MIDDLE
    for i, module in enumerate(bgc_variant.modules):
        module_start_state_idx = module_idx_to_start_state_idx[i]
        next_gene_start_state_idx = module_idx_to_start_state_idx.get(gene_intervals[module.gene_id][1] + 1, -1)
        next_fragment_start_state_idx = module_idx_to_start_state_idx.get(fragment_intervals[module.fragment_idx][1] + 1, None)
        if (i == 0 and next_gene_start_state_idx != -1) or module.gene_id != bgc_variant.modules[i - 1].gene_id:
            adj_list[module_start_state_idx].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.SKIP_GENE,
                                                                    log_prob=0,  # TODO: fill in
                                                                    to=next_gene_start_state_idx))
        if i not in (0, len(bgc_variant.modules) - 1) and module.fragment_idx != bgc_variant.modules[i - 1].fragment_idx:
            adj_list[module_start_state_idx].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.SKIP_FRAGMENT,
                                                                    log_prob=0,  # TODO: fill in
                                                                    to=next_fragment_start_state_idx))




def add_skip_at_the_beginning(states: List[DetailedHMMState],
                              adj_list: Dict[int, List[DetailedHMMEdge]],
                              module_idx_to_start_state_idx: Dict[int, int],
                              fragment_intervals: Dict[int, Tuple[int, int]],
                              initial_state: int = 0) -> None:
    # SKIP AT THE BEGINNING
    num_fragments = len(fragment_intervals)
    fst_fragment_skip_idx = len(states)
    for fragment_idx in range(num_fragments - 1):  # -1 all fragments can't be skipped
        states.append(DetailedHMMState(state_type=DetailedHMMStateType.SKIP_FRAGMENT_AT_BEGINNING,
                                       emissions={}))
        state_idx = len(states) - 1
        if fragment_idx == 0:
            prev_state_idx = initial_state  # initial state
            back_edge_type = DetailedHMMEdgeType.START_SKIPPING_AT_START
        else:
            prev_state_idx = len(states) - 2
            back_edge_type = DetailedHMMEdgeType.CONTINUE_SKIPPING_AT_START
        adj_list[prev_state_idx].append(DetailedHMMEdge(edge_type=back_edge_type,
                                                        log_prob=0,  # TODO: fill in
                                                        to=state_idx))
        next_module_state_idx = module_idx_to_start_state_idx[fragment_intervals[fragment_idx + 1][0]]
        adj_list[state_idx].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.END_SKIPPING_AT_START,
                                                   log_prob=0,  # TODO: fill in
                                                   to=next_module_state_idx))


def add_skip_at_the_end(states: List[DetailedHMMState],
                        adj_list: Dict[int, List[DetailedHMMEdge]],
                        module_idx_to_start_state_idx: Dict[int, int],
                        fragment_intervals: Dict[int, Tuple[int, int]]) -> None:
    # SKIP AT THE END
    num_fragments = len(fragment_intervals)
    for fragment_idx in range(1, num_fragments):  # start from 1 because all fragments can't be skipped
        states.append(DetailedHMMState(state_type=DetailedHMMStateType.SKIP_FRAGMENT_END,
                                       emissions={}))
        state_idx = len(states) - 1
        prev_module_state_idx = module_idx_to_start_state_idx[fragment_intervals[fragment_idx][0]]
        adj_list[prev_module_state_idx].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.START_SKIPPING_AT_END,
                                                               log_prob=0,  # TODO: fill in
                                                               to=state_idx))
        next_state_idx = state_idx + 1  # the next for the last state is naturally the final state, but that looks like a hack, maybe need to redo
        adj_list[state_idx].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.CONTINUE_SKIPPING_AT_END,
                                                   log_prob=0,  # TODO: fill in
                                                   to=next_state_idx))


# TODO: module and gene iterations!
def bgc_variant_to_detailed_hmm(cls,
                                bgc_variant: BGC_Variant):
                                #scoring_helper: ScoringHelper):  # -> DetailedHMM:
    gene_intervals = get_genes_intervals(bgc_variant.modules)
    fragment_intervals = get_fragments_intervals(bgc_variant.modules)

    states = [DetailedHMMState(state_type=DetailedHMMStateType.INITIAL, emissions={})]
    adj_list = defaultdict(list)
    start_state_idx = 0
    final_state_idx = -1  # a bit hacky, but it's fine

    # ADD MODULES' STATES
    module_idx_to_start_state_idx = add_modules(states, adj_list, bgc_variant, final_state_idx)
    state_idx_to_module_idx = {v: k for k, v in module_idx_to_start_state_idx.items()}
    adj_list[start_state_idx].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.START_MATCHING,
                                                    log_prob=0,  # TODO: fill in
                                                    to=module_idx_to_start_state_idx[0]))

    # ADD SKIP GENE AND FRAGMENTS IN THE MIDDLE
    add_skip_gene_and_fragments_middle(adj_list,
                                       bgc_variant,
                                       gene_intervals,
                                       fragment_intervals,
                                       module_idx_to_start_state_idx)

    # SKIP AT THE BEGINNING
    add_skip_at_the_beginning(states, adj_list, module_idx_to_start_state_idx, fragment_intervals)

    # SKIP AT THE END
    add_skip_at_the_end(states, adj_list, module_idx_to_start_state_idx, fragment_intervals)

    # ADD FINAL STATE
    states.append(DetailedHMMState(state_type=DetailedHMMStateType.FINAL, emissions={}))
    final_state_idx = len(states) - 1
    # substitute -1 with final state index
    for state_idx, edges in adj_list.items():
        for i in range(len(edges)):
            if edges[i].to == -1:
                edges[i] = edges[i]._replace(to=final_state_idx)

    return cls(states=states,
               adj_list=adj_list,
               bgc_variant=bgc_variant,
               state_idx_to_module_idx=state_idx_to_module_idx)
