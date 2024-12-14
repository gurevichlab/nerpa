from __future__ import annotations
from typing import List, Tuple, NamedTuple, Dict, Optional
from dataclasses import dataclass
from enum import Enum, auto

from src.data_types import NRP_Monomer, BGC_Module, BGC_Variant, GeneId
from src.matching.matcher_viterbi_types import DetailedHMMEdgeType, DetailedHMMStateType, DetailedHMMState, DetailedHMMEdge, HMM
from src.matching.auxilary import get_genes_intervals, get_fragments_intervals, get_emissions
from collections import defaultdict
from itertools import pairwise
from functools import partial


def gene_len(gene_id: GeneId, genes_intervals: Dict[GeneId, Tuple[int, int]]) -> int:
    return genes_intervals[gene_id][1] - genes_intervals[gene_id][0] + 1


def num_genes_in_fragment(fragment_idx: int, fragment_intervals: Dict[int, Tuple[int, int]],
                          bgc_variant: BGC_Variant) -> int:
    return len(set(bgc_variant.modules[i].gene_id for i in range(fragment_intervals[fragment_idx][0],
                                                                 fragment_intervals[fragment_idx][1] + 1)))


def add_modules(states: List[DetailedHMMState],
                adj_list: Dict[int, List[DetailedHMMEdge]],
                bgc_variant: BGC_Variant,
                gene_intervals: Dict[GeneId, Tuple[int, int]],
                initial_state_idx: int = 0,
                final_state_idx: int = -1) -> Dict[int, int]:  # module_idx -> module_start_state_idx
    _gene_len = partial(gene_len, genes_intervals=gene_intervals)
    module_idx_to_start_state_idx = {i: 3 * i + 1 for i in range(len(bgc_variant.modules))}

    # MODULE STATES: MODULE_START -> MATCH -> INSERT_MONOMER -> MODULE_START -> ...
    for i, module in enumerate(bgc_variant.modules):
        # add states
        states.append(DetailedHMMState(state_type=DetailedHMMStateType.MODULE_START,
                                       emissions={}))
        states.append(DetailedHMMState(state_type=DetailedHMMStateType.MATCH,
                                       emissions=get_emissions(module)))
        states.append(DetailedHMMState(state_type=DetailedHMMStateType.INSERT,
                                       emissions={}))  # TODO: fill in

        # add edges
        start_idx, match_idx, insert_idx = 3 * i + 1, 3 * i + 2, 3 * i + 3
        next_module_idx = 3 * i + 4 if i + 1 < len(bgc_variant.modules) else final_state_idx
        prev_module_idx = 3 * i - 2 if i - 1 >= 0 else initial_state_idx
        adj_list[start_idx].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.MATCH,
                                                   log_prob=0,  # TODO: fill in
                                                   to=match_idx))
        if _gene_len(module.gene_id) > 1:
            adj_list[start_idx].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.SKIP_MODULE,
                                                       log_prob=0,  # TODO: fill in
                                                       to=next_module_idx))
        adj_list[match_idx] = [DetailedHMMEdge(edge_type=DetailedHMMEdgeType.NO_INSERTIONS,
                                               log_prob=0,  # TODO: fill in
                                               to=next_module_idx),
                               DetailedHMMEdge(edge_type=DetailedHMMEdgeType.START_INSERTING,
                                               log_prob=0,  # TODO: fill in
                                               to=insert_idx)]
        adj_list[insert_idx] = [DetailedHMMEdge(edge_type=DetailedHMMEdgeType.INSERT,
                                                log_prob=0,  # TODO: fill in
                                                to=insert_idx),
                                DetailedHMMEdge(edge_type=DetailedHMMEdgeType.END_INSERTING,
                                                log_prob=0,  # TODO: fill in
                                                to=next_module_idx)]
        if module.iterative_module:
            adj_list[next_module_idx].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.ITERATE_MODULE,
                                                             log_prob=0,  # TODO: fill in
                                                             to=start_idx))
        if gene_intervals[module.gene_id][1] == i:  # if the last module of the gene
            fst_module_in_gene_idx = gene_intervals[module.gene_id][0]
            adj_list[next_module_idx].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.ITERATE_GENE,
                                                             log_prob=0,  # TODO: fill in
                                                             to=module_idx_to_start_state_idx[fst_module_in_gene_idx]))
    return module_idx_to_start_state_idx


def add_skip_gene_and_fragments_middle(adj_list: Dict[int, List[DetailedHMMEdge]],
                                       bgc_variant: BGC_Variant,
                                       gene_intervals: Dict[GeneId, Tuple[int, int]],
                                       fragment_intervals: Dict[int, Tuple[int, int]],
                                       module_idx_to_start_state_idx: Dict[int, int]) -> None:
    # SKIP GENES AND FRAGMENTS IN THE MIDDLE
    for i, module in enumerate(bgc_variant.modules[1:], start=1):  # skipping the first gene/fragment is handled separately
        module_start_state_idx = module_idx_to_start_state_idx[i]
        next_gene_start_state_idx = module_idx_to_start_state_idx.get(gene_intervals[module.gene_id][1] + 1, -1)
        next_fragment_start_state_idx = module_idx_to_start_state_idx.get(fragment_intervals[module.fragment_idx][1] + 1, None)
        if module != bgc_variant.modules[i - 1].gene_id:  # module is the first in the gene
            adj_list[module_start_state_idx].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.SKIP_GENE,
                                                                    log_prob=0,  # TODO: fill in
                                                                    to=next_gene_start_state_idx))
        if (module.fragment_idx != bgc_variant.modules[i - 1].fragment_idx   # module is the first in the fragment
                and next_fragment_start_state_idx is not None):  # skipping the last fragment is handled separately
            adj_list[module_start_state_idx].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.SKIP_FRAGMENT,
                                                                    log_prob=0,  # TODO: fill in
                                                                    to=next_fragment_start_state_idx))


def add_skip_at_the_beginning(states: List[DetailedHMMState],
                              adj_list: Dict[int, List[DetailedHMMEdge]],
                              bgc_variant: BGC_Variant,
                              module_idx_to_start_state_idx: Dict[int, int],
                              genes_intervals: Dict[GeneId, Tuple[int, int]],
                              fragment_intervals: Dict[int, Tuple[int, int]],
                              initial_state: int = 0) -> None:
    _gene_len = partial(gene_len, genes_intervals=genes_intervals)
    _num_genes_in_fragment = partial(num_genes_in_fragment, fragment_intervals=fragment_intervals, bgc_variant=bgc_variant)

    # 1. Add inserts before modules
    module_idx_to_insert_state_idx = {}
    for i, module in enumerate(bgc_variant.modules):
        states.append(DetailedHMMState(state_type=DetailedHMMStateType.INSERT_AT_START,
                                       emissions={}))  # TODO: fill emissions
        module_idx_to_insert_state_idx[i] = len(states) - 1
        adj_list[len(states) - 1].extend([DetailedHMMEdge(edge_type=DetailedHMMEdgeType.START_MATCHING,
                                                         log_prob=0,  # TODO: fill in
                                                         to=module_idx_to_start_state_idx[i]),
                                         DetailedHMMEdge(edge_type=DetailedHMMEdgeType.INSERT,
                                                         log_prob=0,  # TODO: fill in
                                                         to=len(states) - 1)])

    # 2. No skips at the beginning
    adj_list[initial_state].extend([DetailedHMMEdge(edge_type=DetailedHMMEdgeType.START_INSERTING_AT_START,
                                                    log_prob=0,  # TODO: fill in
                                                    to=module_idx_to_insert_state_idx[0]),
                                    DetailedHMMEdge(edge_type=DetailedHMMEdgeType.START_MATCHING,
                                                    log_prob=0,  # TODO: fill in
                                                    to=module_idx_to_start_state_idx[0])])

    # 3. Skip modules at the beginning
    module_idx_to_module_skip_state_idx = {}
    for i, module in enumerate(bgc_variant.modules):  # skipping fst_module
        if i == len(bgc_variant.modules) - 1 or bgc_variant.modules[i].gene_id != bgc_variant.modules[i + 1].gene_id:
            continue  # fst_module is the last module of the gene -- skip the gene instead

        states.append(DetailedHMMState(state_type=DetailedHMMStateType.SKIP_MODULE_AT_START,
                                       emissions={}))
        module_idx_to_module_skip_state_idx[i] = len(states) - 1
        adj_list[len(states) - 1].extend([DetailedHMMEdge(edge_type=DetailedHMMEdgeType.START_MATCHING,
                                                          log_prob=0,  # TODO: fill in
                                                          to=module_idx_to_start_state_idx[i + 1]),
                                          DetailedHMMEdge(edge_type=DetailedHMMEdgeType.START_INSERTING_AT_START,
                                                          log_prob=0,  # TODO: fill in
                                                          to=module_idx_to_insert_state_idx[i + 1])])
        if i > 0 and bgc_variant.modules[i].gene_id == bgc_variant.modules[i - 1].gene_id:
            prev_module_skip_state_idx = module_idx_to_module_skip_state_idx[i - 1]
            adj_list[prev_module_skip_state_idx].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.SKIP_MODULE,
                                                                        log_prob=0,  # TODO: fill in
                                                                        to=len(states) - 1))  # skip the previous module

    if _gene_len(bgc_variant.modules[0].gene_id) > 1:
        adj_list[initial_state].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.START_SKIP_MODULES_AT_START,
                                                       log_prob=0,  # TODO: fill in
                                                       to=module_idx_to_module_skip_state_idx[0]))

    # 4. Skip genes at the beginning
    gene_id_to_gene_skip_state_idx = {}
    for fst_module_idx, last_module_idx in sorted(genes_intervals.values()):
        if (last_module_idx == len(bgc_variant.modules) - 1 or
                bgc_variant.modules[last_module_idx].fragment_idx != bgc_variant.modules[last_module_idx + 1].fragment_idx):
            continue  # gene is the last in the fragment -- skip the fragment instead

        states.append(DetailedHMMState(state_type=DetailedHMMStateType.SKIP_GENE_AT_START,
                                        emissions={}))
        gene_id = bgc_variant.modules[fst_module_idx].gene_id
        gene_id_to_gene_skip_state_idx[gene_id] = len(states) - 1
        adj_list[len(states) - 1].extend([DetailedHMMEdge(edge_type=DetailedHMMEdgeType.START_MATCHING,
                                                          log_prob=0,  # TODO: fill in
                                                          to=module_idx_to_start_state_idx[last_module_idx + 1]),
                                          DetailedHMMEdge(edge_type=DetailedHMMEdgeType.START_INSERTING_AT_START,
                                                          log_prob=0,  # TODO: fill in
                                                          to=module_idx_to_insert_state_idx[last_module_idx + 1])])
        if _gene_len(bgc_variant.modules[last_module_idx + 1].gene_id) > 1:  # if the next gene has one module, skip gene instead
            adj_list[len(states) - 1].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.SKIP_MODULE,
                                                             log_prob=0,
                                                             to=module_idx_to_module_skip_state_idx[last_module_idx + 1]))
        if fst_module_idx > 0 and bgc_variant.modules[fst_module_idx].fragment_idx == bgc_variant.modules[fst_module_idx - 1].fragment_idx:
            prev_gene_id = bgc_variant.modules[fst_module_idx - 1].gene_id
            prev_gene_skip_state_idx = gene_id_to_gene_skip_state_idx[prev_gene_id]
            adj_list[prev_gene_skip_state_idx].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.SKIP_GENE,
                                                                     log_prob=0,  # TODO: fill in
                                                                     to=len(states) - 1))

    if _num_genes_in_fragment(0) > 1:
        adj_list[initial_state].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.START_SKIP_GENES_AT_START,
                                                       log_prob=0,  # TODO: fill in
                                                       to=gene_id_to_gene_skip_state_idx[bgc_variant.modules[0].gene_id]))

    # 5. Skip fragments at the beginning
    fragment_idx_to_fragment_skip_state_idx = {}
    for fst_module_idx, last_module_idx in sorted(fragment_intervals.values())[:-1]:  # last fragment can't be skipped
        states.append(DetailedHMMState(state_type=DetailedHMMStateType.SKIP_FRAGMENT_AT_START,
                                        emissions={}))
        fragment_idx = bgc_variant.modules[fst_module_idx].fragment_idx
        fragment_idx_to_fragment_skip_state_idx[fragment_idx] = len(states) - 1
        adj_list[len(states) - 1].extend([DetailedHMMEdge(edge_type=DetailedHMMEdgeType.START_MATCHING,
                                                          log_prob=0,  # TODO: fill in
                                                          to=module_idx_to_start_state_idx[last_module_idx + 1]),
                                            DetailedHMMEdge(edge_type=DetailedHMMEdgeType.START_INSERTING_AT_START,
                                                            log_prob=0,  # TODO: fill in
                                                            to=module_idx_to_insert_state_idx[last_module_idx + 1])])
        next_gene_id = bgc_variant.modules[last_module_idx + 1].gene_id
        if _num_genes_in_fragment(fragment_idx + 1) > 1:
            adj_list[len(states) - 1].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.SKIP_GENE,
                                                             log_prob=0,  # TODO: fill in
                                                             to=gene_id_to_gene_skip_state_idx[next_gene_id]))
        if _gene_len(next_gene_id) > 1:
            adj_list[len(states) - 1].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.SKIP_MODULE,
                                                             log_prob=0,  # TODO: fill in
                                                             to=module_idx_to_module_skip_state_idx[last_module_idx + 1]))
        if fst_module_idx > 0:
            prev_fragment_idx = bgc_variant.modules[fst_module_idx - 1].fragment_idx
            prev_fragment_skip_state_idx = fragment_idx_to_fragment_skip_state_idx[prev_fragment_idx]
            adj_list[prev_fragment_skip_state_idx].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.SKIP_FRAGMENT,
                                                                        log_prob=0,  # TODO: fill in
                                                                        to=len(states) - 1))

    if len(fragment_intervals) > 1:
        adj_list[initial_state].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.START_SKIP_FRAGMENTS_AT_START,
                                                       log_prob=0,  # TODO: fill in
                                                       to=fragment_idx_to_fragment_skip_state_idx[0]))


def add_skip_at_the_end(states: List[DetailedHMMState],
                        adj_list: Dict[int, List[DetailedHMMEdge]],
                        module_idx_to_start_state_idx: Dict[int, int],
                        fragment_intervals: Dict[int, Tuple[int, int]],
                        final_state_idx = -1) -> None:
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


# TODO:
#  1. module and gene iterations!  # DONE, though right now they are unlimited
#  2. insertions at start
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
    module_idx_to_state_idx = add_modules(states,
                                          adj_list,
                                          bgc_variant,
                                          gene_intervals,
                                          start_state_idx,
                                          final_state_idx)  # module_idx -> module_start_state_idx
    state_idx_to_module_idx = {v: k for k, v in module_idx_to_state_idx.items()}
    adj_list[start_state_idx].append(DetailedHMMEdge(edge_type=DetailedHMMEdgeType.START_MATCHING,
                                                    log_prob=0,  # TODO: fill in
                                                    to=module_idx_to_state_idx[0]))

    # ADD SKIP GENE AND FRAGMENTS IN THE MIDDLE
    add_skip_gene_and_fragments_middle(adj_list,
                                       bgc_variant,
                                       gene_intervals,
                                       fragment_intervals,
                                       module_idx_to_state_idx)

    # SKIP AT THE BEGINNING
    add_skip_at_the_beginning(states,
                              adj_list,
                              bgc_variant,
                              module_idx_to_state_idx,
                              gene_intervals,
                              fragment_intervals)

    # SKIP AT THE END
    add_skip_at_the_end(states, adj_list, module_idx_to_state_idx, fragment_intervals)

    # ADD FINAL STATE
    states.append(DetailedHMMState(state_type=DetailedHMMStateType.FINAL, emissions={}))
    final_state_idx = len(states) - 1
    # substitute -1 with final state index
    if -1 in adj_list:
        adj_list[final_state_idx] = adj_list.pop(-1)
    for state_idx, edges in adj_list.items():
        for i in range(len(edges)):
            if edges[i].to == -1:
                edges[i] = edges[i]._replace(to=final_state_idx)

    return cls(states=states,
               adj_list=adj_list,
               bgc_variant=bgc_variant,
               state_idx_to_module_idx=state_idx_to_module_idx,
               _module_idx_to_state_idx=module_idx_to_state_idx,
               start_state_idx=start_state_idx,
               final_state_idx=final_state_idx)
