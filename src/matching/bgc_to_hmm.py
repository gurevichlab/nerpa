from __future__ import annotations
from typing import List, Tuple, NamedTuple, Dict, Optional
from dataclasses import dataclass
from enum import Enum, auto

from src.antismash_parsing.genomic_context import (
    ModuleGenomicContext,
    GeneGenomicContext,
    FragmentGenomicContext,
    module_features_to_gene_features,
    module_features_to_fragment_features
)
from src.data_types import BGC_Variant, GeneId
from src.matching.hmm_auxiliary_types import DetailedHMMEdgeType, DetailedHMMStateType, DetailedHMMState, DetailedHMMEdge, HMM
from src.matching.genes_fragments_intervals import get_genes_intervals, get_fragments_intervals
from collections import defaultdict
from functools import partial
from collections import OrderedDict
from src.matching.hmm_edge_weights import get_edge_weights
from src.matching.hmm_scoring_helper import HMMHelper
import networkx as nx

def make_edge(edge_type: DetailedHMMEdgeType,
              log_prob: float,
              bgc_variant: BGC_Variant,
              fst_module_idx: Optional[int] = None) -> DetailedHMMEdge:
    if fst_module_idx is None:
        return DetailedHMMEdge(edge_type=edge_type,
                               weight=log_prob,
                               genomic_context=None,
                               edge_key=(bgc_variant.genome_id, bgc_variant.bgc_idx))

    fst_module = bgc_variant.modules[fst_module_idx]

    if edge_type == DetailedHMMEdgeType.SKIP_GENE:
        gene_id = bgc_variant.modules[fst_module_idx].gene_id
        lst_module_idx = next((i for i in range(fst_module_idx, len(bgc_variant.modules) - 1)
                               if bgc_variant.modules[i + 1].gene_id != gene_id),
                              len(bgc_variant.modules) - 1)
        gene_context = module_features_to_gene_features(bgc_variant.modules[fst_module_idx].genomic_context,
                                                        bgc_variant.modules[lst_module_idx].genomic_context)
        return DetailedHMMEdge(edge_type=edge_type,
                               weight=log_prob,
                               genomic_context=gene_context,
                               edge_key=(bgc_variant.genome_id, bgc_variant.bgc_idx,
                                         fst_module.gene_id))
    if edge_type in (
            DetailedHMMEdgeType.SKIP_FRAGMENT,
            DetailedHMMEdgeType.SKIP_FRAGMENT_AT_START,
            DetailedHMMEdgeType.SKIP_FRAGMENT_AT_END):
        fragment_idx = bgc_variant.modules[fst_module_idx].fragment_idx
        lst_module_idx = next((i for i in range(fst_module_idx, len(bgc_variant.modules) - 1)
                               if bgc_variant.modules[i + 1].fragment_idx != fragment_idx),
                                len(bgc_variant.modules) - 1)
        fragment_context = module_features_to_fragment_features(bgc_variant.modules[fst_module_idx].genomic_context,
                                                                bgc_variant.modules[lst_module_idx].genomic_context)
        genes_ids = list(OrderedDict.fromkeys(bgc_variant.modules[i].gene_id
                                              for i in range(fst_module_idx, lst_module_idx + 1)))
        return DetailedHMMEdge(edge_type=edge_type,
                               weight=log_prob,
                               genomic_context=fragment_context,
                               edge_key=(bgc_variant.genome_id, bgc_variant.bgc_idx,
                                         *genes_ids))
    module_context = bgc_variant.modules[fst_module_idx].genomic_context
    return DetailedHMMEdge(edge_type=edge_type,
                           weight=log_prob,
                           genomic_context=module_context,
                           edge_key=(bgc_variant.genome_id, bgc_variant.bgc_idx,
                                     fst_module.gene_id, fst_module.a_domain_idx))


def gene_len(gene_id: GeneId, genes_intervals: Dict[GeneId, Tuple[int, int]]) -> int:
    return genes_intervals[gene_id][1] - genes_intervals[gene_id][0] + 1


def num_genes_in_fragment(fragment_idx: int, fragment_intervals: Dict[int, Tuple[int, int]],
                          bgc_variant: BGC_Variant) -> int:
    return len(set(bgc_variant.modules[i].gene_id for i in range(fragment_intervals[fragment_idx][0],
                                                                 fragment_intervals[fragment_idx][1] + 1)))


def add_modules(states: List[DetailedHMMState],
                adj_list: Dict[int, Dict[int, DetailedHMMEdge]],
                bgc_variant: BGC_Variant,
                hmm_helper: HMMHelper,
                gene_intervals: Dict[GeneId, Tuple[int, int]],
                initial_state_idx: int = 0,
                final_state_idx: int = -1) -> Dict[int, int]:  # module_idx -> module_start_state_idx
    _gene_len = partial(gene_len, genes_intervals=gene_intervals)
    _make_edge = partial(make_edge, bgc_variant=bgc_variant)
    _get_emissions = partial(hmm_helper.get_emissions, pks_domains_in_bgc=bgc_variant.has_pks_domains())
    _get_insert_emissions = partial(hmm_helper.get_insert_emissions, pks_domains_in_bgc=bgc_variant.has_pks_domains())

    module_idx_to_start_state_idx = {i: 3 * i + 1 for i in range(len(bgc_variant.modules))}

    # MODULE STATES: MODULE_START -> MATCH -> INSERT_MONOMER -> MODULE_START -> ...
    for i, module in enumerate(bgc_variant.modules):
        # add states
        states.extend([
            DetailedHMMState(state_type=DetailedHMMStateType.MODULE_START,
                             emissions={}),
            DetailedHMMState(state_type=DetailedHMMStateType.MATCH,
                             emissions=_get_emissions(module)),
            DetailedHMMState(state_type=DetailedHMMStateType.INSERT,
                             emissions=_get_insert_emissions(module)),
        ])

        # add edges
        start_idx, match_idx, insert_idx = 3 * i + 1, 3 * i + 2, 3 * i + 3
        next_module_idx = 3 * i + 4 if i + 1 < len(bgc_variant.modules) else final_state_idx
        adj_list[start_idx][match_idx] = _make_edge(edge_type=DetailedHMMEdgeType.MATCH,
                                                    log_prob=0,  # TODO: fill in
                                                    fst_module_idx=i)

        if _gene_len(module.gene_id) > 1:
            adj_list[start_idx][next_module_idx] = _make_edge(edge_type=DetailedHMMEdgeType.SKIP_MODULE,
                                                              log_prob=0,  # TODO: fill in
                                                              fst_module_idx=i)
        adj_list[match_idx][next_module_idx] = _make_edge(edge_type=DetailedHMMEdgeType.NO_INSERTIONS,
                                                          log_prob=0,  # TODO: fill in
                                                          fst_module_idx=i)
        adj_list[match_idx][insert_idx] = _make_edge(edge_type=DetailedHMMEdgeType.START_INSERTING,
                                                     log_prob=0,  # TODO: fill in
                                                     fst_module_idx=i)
        adj_list[insert_idx][insert_idx] = _make_edge(edge_type=DetailedHMMEdgeType.INSERT,
                                                      log_prob=0,  # TODO: fill in
                                                      fst_module_idx=None)
        adj_list[insert_idx][next_module_idx] = _make_edge(edge_type=DetailedHMMEdgeType.END_INSERTING,
                                                           log_prob=0,  # TODO: fill in
                                                           fst_module_idx=None)
    return module_idx_to_start_state_idx


def add_modules_and_genes_iterations(adj_list: Dict[int, Dict[int, DetailedHMMEdge]],
                                     bgc_variant: BGC_Variant,
                                     gene_intervals: Dict[GeneId, Tuple[int, int]],
                                     module_idx_to_state_idx: Dict[int, int],
                                     final_state_idx: int) -> None:
    _make_edge = partial(make_edge, bgc_variant=bgc_variant)

    for module_idx, module in enumerate(bgc_variant.modules):
        start_idx = module_idx_to_state_idx[module_idx]
        next_module_start_idx = module_idx_to_state_idx[module_idx + 1] \
            if module_idx + 1 < len(bgc_variant.modules) else final_state_idx
        if module.iterative_module:
            adj_list[next_module_start_idx][start_idx] = _make_edge(edge_type=DetailedHMMEdgeType.ITERATE_MODULE,
                                                                    log_prob=0,  # TODO: fill in
                                                                    fst_module_idx=None)
        if module.iterative_gene and gene_intervals[module.gene_id][1] == module_idx:  # if the last module of the gene
            fst_module_in_gene_idx = gene_intervals[module.gene_id][0]
            edge_end = module_idx_to_state_idx[fst_module_in_gene_idx]
            adj_list[next_module_start_idx][edge_end] = _make_edge(edge_type=DetailedHMMEdgeType.ITERATE_GENE,
                                                                   log_prob=0,  # TODO: fill in
                                                                   fst_module_idx=None)


def add_skip_gene_and_fragments_middle(adj_list: Dict[int, Dict[int, DetailedHMMEdge]],
                                       bgc_variant: BGC_Variant,
                                       gene_intervals: Dict[GeneId, Tuple[int, int]],
                                       fragment_intervals: Dict[int, Tuple[int, int]],
                                       module_idx_to_start_state_idx: Dict[int, int],
                                       final_state_idx: int) -> None:
    _make_edge = partial(make_edge, bgc_variant=bgc_variant)

    # SKIP GENES AND FRAGMENTS IN THE MIDDLE
    for i, module in enumerate(bgc_variant.modules[1:], start=1):  # skipping the first gene/fragment is handled separately
        module_start_state_idx = module_idx_to_start_state_idx[i]
        next_gene_start_state_idx = module_idx_to_start_state_idx.get(gene_intervals[module.gene_id][1] + 1, final_state_idx)
        next_fragment_start_state_idx = module_idx_to_start_state_idx.get(fragment_intervals[module.fragment_idx][1] + 1, None)
        if module.gene_id != bgc_variant.modules[i - 1].gene_id:  # module is the first in the gene
            adj_list[module_start_state_idx][next_gene_start_state_idx] = _make_edge(edge_type=DetailedHMMEdgeType.SKIP_GENE,
                                                                                     log_prob=0,  # TODO: fill in
                                                                                     fst_module_idx=i)
        if (module.fragment_idx != bgc_variant.modules[i - 1].fragment_idx   # module is the first in the fragment
                and next_fragment_start_state_idx is not None):  # skipping the last fragment is handled separately
            adj_list[module_start_state_idx][next_fragment_start_state_idx] = _make_edge(edge_type=DetailedHMMEdgeType.SKIP_FRAGMENT,
                                                                                         log_prob=0,  # TODO: fill in
                                                                                         fst_module_idx=i)


def add_skip_at_the_beginning(states: List[DetailedHMMState],
                              adj_list: Dict[int, Dict[int, DetailedHMMEdge]],
                              bgc_variant: BGC_Variant,
                              hmm_helper: HMMHelper,
                              module_idx_to_start_state_idx: Dict[int, int],
                              genes_intervals: Dict[GeneId, Tuple[int, int]],
                              fragment_intervals: Dict[int, Tuple[int, int]],
                              initial_state: int = 0) -> None:
    _make_edge = partial(make_edge, bgc_variant=bgc_variant)
    _gene_len = partial(gene_len, genes_intervals=genes_intervals)
    _num_genes_in_fragment = partial(num_genes_in_fragment, fragment_intervals=fragment_intervals, bgc_variant=bgc_variant)
    _get_insert_at_start_emissions = partial(hmm_helper.get_insert_at_start_emissions, pks_domains_in_bgc=bgc_variant.has_pks_domains())

    # 1. Add inserts before modules
    module_idx_to_insert_state_idx = {}
    for i, module in enumerate(bgc_variant.modules):
        states.append(DetailedHMMState(state_type=DetailedHMMStateType.INSERT_AT_START,
                                       emissions=_get_insert_at_start_emissions(module)))
        module_idx_to_insert_state_idx[i] = len(states) - 1
        adj_list[len(states) - 1][module_idx_to_start_state_idx[i]] = _make_edge(edge_type=DetailedHMMEdgeType.START_MATCHING,
                                                                                 log_prob=0,  # TODO: fill in
                                                                                 fst_module_idx=None)
        adj_list[len(states) - 1][len(states) - 1] = _make_edge(edge_type=DetailedHMMEdgeType.INSERT,
                                                                log_prob=0,  # TODO: fill in
                                                                fst_module_idx=None)

    # 2. No skips at the beginning
    adj_list[initial_state][module_idx_to_insert_state_idx[0]] = _make_edge(edge_type=DetailedHMMEdgeType.START_INSERTING_AT_START,
                                                                            log_prob=0,  # TODO: fill in
                                                                            fst_module_idx=0)
    adj_list[initial_state][module_idx_to_start_state_idx[0]] = _make_edge(edge_type=DetailedHMMEdgeType.START_MATCHING,
                                                                           log_prob=0,  # TODO: fill in
                                                                           fst_module_idx=0)

    # 3. Skip modules at the beginning
    module_idx_to_module_skip_state_idx = {}
    for i, module in enumerate(bgc_variant.modules):  # skipping modules[:i+1]
        if i == len(bgc_variant.modules) - 1 or bgc_variant.modules[i].gene_id != bgc_variant.modules[i + 1].gene_id:
            continue  # module is the last module of the gene -- skip the gene instead

        states.append(DetailedHMMState(state_type=DetailedHMMStateType.SKIP_MODULE_AT_START,
                                       emissions={}))
        module_idx_to_module_skip_state_idx[i] = len(states) - 1
        adj_list[len(states) - 1][module_idx_to_start_state_idx[i + 1]] = _make_edge(edge_type=DetailedHMMEdgeType.START_MATCHING,
                                                                                     log_prob=0,  # TODO: fill in
                                                                                     fst_module_idx=i + 1)
        adj_list[len(states) - 1][module_idx_to_insert_state_idx[i + 1]] = _make_edge(edge_type=DetailedHMMEdgeType.START_INSERTING_AT_START,
                                                                                      log_prob=0,  # TODO: fill in
                                                                                      fst_module_idx=i + 1)
        if i > 0 and bgc_variant.modules[i].gene_id == bgc_variant.modules[i - 1].gene_id:
            prev_module_skip_state_idx = module_idx_to_module_skip_state_idx[i - 1]
            adj_list[prev_module_skip_state_idx][len(states) - 1] = _make_edge(edge_type=DetailedHMMEdgeType.SKIP_MODULE,
                                                                               log_prob=0,  # TODO: fill in
                                                                               fst_module_idx=i)  # skip the previous module

    if _gene_len(bgc_variant.modules[0].gene_id) > 1:
        adj_list[initial_state][module_idx_to_module_skip_state_idx[0]] = _make_edge(edge_type=DetailedHMMEdgeType.START_SKIP_MODULES_AT_START,
                                                                                     log_prob=0,  # TODO: fill in
                                                                                     fst_module_idx=None)

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
        adj_list[len(states) - 1][module_idx_to_start_state_idx[last_module_idx + 1]] = _make_edge(edge_type=DetailedHMMEdgeType.START_MATCHING,
                                                                                                   log_prob=0,  # TODO: fill in
                                                                                                   fst_module_idx=last_module_idx + 1)
        adj_list[len(states) - 1][module_idx_to_insert_state_idx[last_module_idx + 1]] = _make_edge(edge_type=DetailedHMMEdgeType.START_INSERTING_AT_START,
                                                                                                    log_prob=0,  # TODO: fill in
                                                                                                    fst_module_idx=last_module_idx + 1)
        if _gene_len(bgc_variant.modules[last_module_idx + 1].gene_id) > 1:  # if the next gene has one module, skip gene instead
            adj_list[len(states) - 1][module_idx_to_module_skip_state_idx[last_module_idx + 1]] = _make_edge(edge_type=DetailedHMMEdgeType.SKIP_MODULE,
                                                                                                             log_prob=0,  # TODO: fill in
                                                                                                             fst_module_idx=last_module_idx + 1)
        if fst_module_idx > 0 and bgc_variant.modules[fst_module_idx].fragment_idx == bgc_variant.modules[fst_module_idx - 1].fragment_idx:
            prev_gene_id = bgc_variant.modules[fst_module_idx - 1].gene_id
            prev_gene_skip_state_idx = gene_id_to_gene_skip_state_idx[prev_gene_id]
            adj_list[prev_gene_skip_state_idx][len(states) - 1] = _make_edge(edge_type=DetailedHMMEdgeType.SKIP_GENE,
                                                                             log_prob=0,  # TODO: fill in
                                                                             fst_module_idx=fst_module_idx)

    if _num_genes_in_fragment(0) > 1:
        adj_list[initial_state][gene_id_to_gene_skip_state_idx[bgc_variant.modules[0].gene_id]] = _make_edge(edge_type=DetailedHMMEdgeType.START_SKIP_GENES_AT_START,
                                                                                                             log_prob=0,  # TODO: fill in
                                                                                                             fst_module_idx=0)

    # 5. Skip fragments at the beginning
    fragment_idx_to_fragment_skip_state_idx = {}
    for fst_module_idx, last_module_idx in sorted(fragment_intervals.values())[:-1]:  # last fragment can't be skipped
        states.append(DetailedHMMState(state_type=DetailedHMMStateType.SKIP_FRAGMENT_AT_START,
                                        emissions={}))
        fragment_idx = bgc_variant.modules[fst_module_idx].fragment_idx
        fragment_idx_to_fragment_skip_state_idx[fragment_idx] = len(states) - 1
        adj_list[len(states) - 1][module_idx_to_start_state_idx[last_module_idx + 1]] = _make_edge(edge_type=DetailedHMMEdgeType.START_MATCHING,
                                                                                                    log_prob=0,  # TODO: fill in
                                                                                                   fst_module_idx=last_module_idx + 1)
        adj_list[len(states) - 1][module_idx_to_insert_state_idx[last_module_idx + 1]] = _make_edge(edge_type=DetailedHMMEdgeType.START_INSERTING_AT_START,
                                                                                                    log_prob=0,  # TODO: fill in
                                                                                                    fst_module_idx=last_module_idx + 1)
        next_gene_id = bgc_variant.modules[last_module_idx + 1].gene_id
        if _num_genes_in_fragment(fragment_idx + 1) > 1:
            adj_list[len(states) - 1][gene_id_to_gene_skip_state_idx[next_gene_id]] = _make_edge(edge_type=DetailedHMMEdgeType.SKIP_GENE,
                                                                                                 log_prob=0,  # TODO: fill in
                                                                                                 fst_module_idx=fst_module_idx)
        if _gene_len(next_gene_id) > 1:
            adj_list[len(states) - 1][module_idx_to_module_skip_state_idx[last_module_idx + 1]] = _make_edge(edge_type=DetailedHMMEdgeType.SKIP_MODULE,
                                                                                                             log_prob=0,  # TODO: fill in
                                                                                                             fst_module_idx=last_module_idx + 1)
        if fst_module_idx > 0:
            prev_fragment_idx = bgc_variant.modules[fst_module_idx - 1].fragment_idx
            prev_fragment_skip_state_idx = fragment_idx_to_fragment_skip_state_idx[prev_fragment_idx]
            adj_list[prev_fragment_skip_state_idx][len(states) - 1] = _make_edge(edge_type=DetailedHMMEdgeType.SKIP_FRAGMENT,
                                                                                 log_prob=0,  # TODO: fill in
                                                                                 fst_module_idx=fst_module_idx)  # I'm not sure if this is correct but I am lost at this point

    if len(fragment_intervals) > 1:
        adj_list[initial_state][fragment_idx_to_fragment_skip_state_idx[0]] = _make_edge(edge_type=DetailedHMMEdgeType.START_SKIP_FRAGMENTS_AT_START,
                                                                                         log_prob=0,  # TODO: fill in
                                                                                         fst_module_idx=0)


def add_skip_at_the_end(states: List[DetailedHMMState],
                        adj_list: Dict[int, Dict[int, DetailedHMMEdge]],
                        bgc_variant: BGC_Variant,
                        module_idx_to_start_state_idx: Dict[int, int],
                        fragment_intervals: Dict[int, Tuple[int, int]],
                        final_state_idx = -1) -> None:
    _make_edge = partial(make_edge, bgc_variant=bgc_variant)
    # SKIP AT THE END
    num_fragments = len(fragment_intervals)
    for fragment_idx in range(1, num_fragments):  # start from 1 because all fragments can't be skipped
        states.append(DetailedHMMState(state_type=DetailedHMMStateType.SKIP_FRAGMENT_END,
                                       emissions={}))
        state_idx = len(states) - 1
        prev_module_idx = fragment_intervals[fragment_idx][0]
        prev_module_state_idx = module_idx_to_start_state_idx[prev_module_idx]
        adj_list[prev_module_state_idx][state_idx] = _make_edge(edge_type=DetailedHMMEdgeType.START_SKIPPING_AT_END,
                                                                log_prob=0,  # TODO: fill in
                                                                fst_module_idx=prev_module_idx)
        next_state_idx = state_idx + 1  # the next for the last state is naturally the final state, but that looks like a hack, maybe need to redo
        adj_list[state_idx][next_state_idx] = _make_edge(edge_type=DetailedHMMEdgeType.SKIP_FRAGMENT_AT_END,
                                                         log_prob=0,  # TODO: fill in
                                                         fst_module_idx=None)


def topsort_states(states: List[DetailedHMMState],
                   adj_list: Dict[int, Dict[int, DetailedHMMEdge]],
                   module_idx_to_state_idx: Dict[int, int],
                   state_idx_to_module_idx: Dict[int, Optional[int]]) -> None:
    """Sort states topologically using networkx."""

    # Create a directed graph from the adjacency list
    G = nx.DiGraph()
    for state, neighbors in adj_list.items():
        for neighbor in neighbors:
            if state != neighbor:
                G.add_edge(state, neighbor)

    # q: find cycle in G
    try:
        cycle = nx.find_cycle(G)
        raise ValueError(f"Cycle detected in the HMM state graph: {cycle}")
    except nx.NetworkXNoCycle:
        pass

    sorted_states = list(nx.topological_sort(G))

    # Create a new mapping for state indices based on topological order
    new_state_mapping = {old_idx: new_idx for new_idx, old_idx in enumerate(sorted_states)}
    states[:] = [states[old_idx] for old_idx in sorted_states]

    # Preserve the current state-to-module mapping before clearing
    old_state_idx_to_module_idx = state_idx_to_module_idx.copy()

    # Clear and update module index mappings
    state_idx_to_module_idx.clear()
    module_idx_to_state_idx.clear()

    for new_idx, old_idx in enumerate(sorted_states):
        if old_idx in old_state_idx_to_module_idx:
            module_idx_to_state_idx[old_state_idx_to_module_idx[old_idx]] = new_idx
            state_idx_to_module_idx[new_idx] = old_state_idx_to_module_idx[old_idx]

    # Rebuild adjacency list with new state indices
    new_adj_list = defaultdict(dict)
    for old_idx, neighbors in adj_list.items():
        new_idx = new_state_mapping[old_idx]
        for neighbor, edge in neighbors.items():
            new_adj_list[new_idx][new_state_mapping[neighbor]] = edge

    adj_list.clear()
    adj_list.update(new_adj_list)


def bgc_variant_to_detailed_hmm(cls,
                                bgc_variant: BGC_Variant,
                                hmm_helper: HMMHelper):
    _make_edge = partial(make_edge, bgc_variant=bgc_variant)
    gene_intervals = get_genes_intervals(bgc_variant.modules)
    fragment_intervals = get_fragments_intervals(bgc_variant.modules)

    states = [DetailedHMMState(state_type=DetailedHMMStateType.INITIAL, emissions={})]
    adj_list = defaultdict(dict)
    start_state_idx = 0
    final_state_idx = -1  # a bit hacky, but it's fine

    # ADD MODULES' STATES
    module_idx_to_state_idx = add_modules(states,
                                          adj_list,
                                          bgc_variant,
                                          hmm_helper,
                                          gene_intervals,
                                          start_state_idx,
                                          final_state_idx)  # module_idx -> module_start_state_idx
    state_idx_to_module_idx = {v: k for k, v in module_idx_to_state_idx.items()}
    adj_list[start_state_idx][module_idx_to_state_idx[0]] = _make_edge(edge_type=DetailedHMMEdgeType.START_MATCHING,
                                                                       log_prob=0,  # TODO: fill in
                                                                       fst_module_idx=0)


    # SKIP AT THE BEGINNING
    add_skip_at_the_beginning(states,
                              adj_list,
                              bgc_variant,
                              hmm_helper,
                              module_idx_to_state_idx,
                              gene_intervals,
                              fragment_intervals)

    # SKIP AT THE END
    add_skip_at_the_end(states,
                        adj_list,
                        bgc_variant,
                        module_idx_to_state_idx,
                        fragment_intervals)

    # ADD FINAL STATE
    states.append(DetailedHMMState(state_type=DetailedHMMStateType.FINAL, emissions={}))
    final_state_idx = len(states) - 1
    # substitute -1 with final state index
    if -1 in adj_list:
        adj_list[final_state_idx] = adj_list.pop(-1)
    for state_idx, neighbors in adj_list.items():
        if -1 in neighbors:
            neighbors[final_state_idx] = neighbors[-1]
            del neighbors[-1]

    # Enumerate states in topological order
    topsort_states(states,
                   adj_list,
                   module_idx_to_state_idx,
                   state_idx_to_module_idx)

    # ADD SKIP GENE AND FRAGMENTS IN THE MIDDLE
    add_skip_gene_and_fragments_middle(adj_list,
                                       bgc_variant,
                                       gene_intervals,
                                       fragment_intervals,
                                       module_idx_to_state_idx,
                                       final_state_idx)

    add_modules_and_genes_iterations(adj_list,
                                     bgc_variant,
                                     gene_intervals,
                                     module_idx_to_state_idx,
                                     final_state_idx)

    # Calculate module_idx -> module_match_state_idx
    match_states = [state_idx for state_idx, state in enumerate(states)
                    if state.state_type == DetailedHMMStateType.MATCH]

    hmm = cls(states=states,
              transitions=adj_list,
              bgc_variant=bgc_variant,
              state_idx_to_module_idx=state_idx_to_module_idx,
              _module_idx_to_state_idx=module_idx_to_state_idx,
              _module_idx_to_match_state_idx=match_states,
              start_state_idx=start_state_idx,
              final_state_idx=final_state_idx,
              hmm_helper=hmm_helper)

    # set edge weights
    edge_weights = hmm_helper.get_edge_weights(hmm)
    for u, v in edge_weights:
       hmm.transitions[u][v] = hmm.transitions[u][v]._replace(weight=edge_weights[(u, v)])

    return hmm
