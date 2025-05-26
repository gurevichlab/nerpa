from __future__ import annotations
from typing import List, Tuple, NamedTuple, Dict, Optional

from src.antismash_parsing.genomic_context import ModuleGenomicContextFeature
from src.data_types import BGC_Variant, GeneId
from src.matching.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    DetailedHMMStateType,
    DetailedHMMState,
    DetailedHMMEdge,
    ModuleLevelEdgeKey,
    GeneLevelEdgeKey,
)
from src.matching.genes_fragments_intervals import get_genes_intervals, get_fragments_intervals
from collections import defaultdict
from functools import partial
from src.matching.hmm_scoring_helper import HMMHelper
import networkx as nx

def make_edge(edge_type: DetailedHMMEdgeType,
              log_prob: float,
              bgc_variant: BGC_Variant,
              fst_module_idx: Optional[int] = None) -> DetailedHMMEdge:
    if fst_module_idx is None:
        return DetailedHMMEdge(edge_type=edge_type,
                               weight=log_prob,
                               genomic_context=(),
                               edge_key=None)

    fst_module = bgc_variant.modules[fst_module_idx]

    if edge_type == DetailedHMMEdgeType.SKIP_GENE:
        return DetailedHMMEdge(edge_type=edge_type,
                               weight=log_prob,
                               genomic_context=(),
                               edge_key=GeneLevelEdgeKey(fst_module.gene_id))

    module_context = bgc_variant.modules[fst_module_idx].genomic_context
    return DetailedHMMEdge(edge_type=edge_type,
                           weight=log_prob,
                           genomic_context=module_context,
                           edge_key=ModuleLevelEdgeKey((fst_module.gene_id, fst_module.a_domain_idx)))


def gene_len(gene_id: GeneId, genes_intervals: Dict[GeneId, Tuple[int, int]]) -> int:
    return genes_intervals[gene_id][1] - genes_intervals[gene_id][0] + 1


def add_modules(states: List[DetailedHMMState],
                adj_list: Dict[int, Dict[int, DetailedHMMEdge]],
                bgc_variant: BGC_Variant,
                hmm_helper: HMMHelper,
                gene_intervals: Dict[GeneId, Tuple[int, int]],
                initial_state_idx: int = 0,
                final_state_idx: int = -1) -> Dict[int, int]:  # module_idx -> module_start_state_idx
    ST = DetailedHMMStateType
    ET = DetailedHMMEdgeType

    _gene_len = partial(gene_len, genes_intervals=gene_intervals)
    _make_edge = partial(make_edge,
                         bgc_variant=bgc_variant,
                         log_prob=None,
                         fst_module_idx=None)
    _get_emissions = partial(hmm_helper.get_emissions,
                             pks_domains_in_bgc=bgc_variant.has_pks_domains())
    _get_insert_emissions = partial(hmm_helper.get_insert_emissions,
                                    pks_domains_in_bgc=bgc_variant.has_pks_domains())

    module_idx_to_start_state_idx = {}

    # MODULE STATES: MODULE_START -> MATCH -> INSERT_MONOMER -> MODULE_START -> ...
    for i, module in enumerate(bgc_variant.modules):
        # add states
        states.extend([
            DetailedHMMState(state_type=ST.MODULE_START, emissions={}),
            DetailedHMMState(state_type=ST.MATCH, emissions=_get_emissions(module)),
            DetailedHMMState(state_type=(ST.INSERT if i < len(bgc_variant.modules) - 1
                                         else ST.INSERT_AT_END),
                             emissions=_get_insert_emissions(module))
        ])

        start_idx, match_idx, insert_idx = range(len(states) - 3, len(states))
        module_idx_to_start_state_idx[i] = start_idx

        iteration = False
        iteration_idx = None
        if module.iterative_module or (module.iterative_gene and gene_intervals[module.gene_id][1] == i):
            iteration = True
            states.append(DetailedHMMState(state_type=ST.CHOOSE_IF_ITERATE, emissions={}))
            iteration_idx = len(states) - 1

        # add edges
        next_module_idx = len(states) if i + 1 < len(bgc_variant.modules) else final_state_idx
        adj_list[start_idx][match_idx] = _make_edge(edge_type=ET.MATCH,
                                                    fst_module_idx=i)
        adj_list[match_idx][insert_idx] = _make_edge(edge_type=(ET.START_INSERTING
                                                                if i < len(bgc_variant.modules) - 1
                                                                else ET.START_INSERTING_AT_END),
                                                     fst_module_idx=i)
        adj_list[insert_idx][insert_idx] = _make_edge(edge_type=(ET.INSERT
                                                                 if i < len(bgc_variant.modules) - 1
                                                                 else ET.INSERT_AT_END),
                                                      fst_module_idx=None)

        #if _gene_len(module.gene_id) > 1:
        if iteration:
            adj_list[start_idx][iteration_idx] = _make_edge(edge_type=ET.SKIP_MODULE,
                                                            fst_module_idx=i)
            adj_list[match_idx][iteration_idx] = _make_edge(edge_type=ET.NO_INSERTIONS,
                                                            fst_module_idx=i)
            adj_list[insert_idx][iteration_idx] = _make_edge(edge_type=ET.END_INSERTING)
            adj_list[iteration_idx][next_module_idx] = _make_edge(edge_type=ET.NO_ITERATION,
                                                                  fst_module_idx=i)
        else:
            adj_list[start_idx][next_module_idx] = _make_edge(edge_type=ET.SKIP_MODULE,
                                                              fst_module_idx=i)
            adj_list[match_idx][next_module_idx] = _make_edge(edge_type=ET.NO_INSERTIONS,
                                                              fst_module_idx=i)
            adj_list[insert_idx][next_module_idx] = _make_edge(edge_type=ET.END_INSERTING)

    return module_idx_to_start_state_idx


def add_modules_and_genes_iterations(states: List[DetailedHMMState],
                                     adj_list: Dict[int, Dict[int, DetailedHMMEdge]],
                                     bgc_variant: BGC_Variant,
                                     gene_intervals: Dict[GeneId, Tuple[int, int]],
                                     module_idx_to_state_idx: Dict[int, int],
                                     final_state_idx: int) -> None:
    _make_edge = partial(make_edge, bgc_variant=bgc_variant)

    for module_idx, module in enumerate(bgc_variant.modules):
        if not module.iterative_module and not (module.iterative_gene and gene_intervals[module.gene_id][1] == module_idx):
            continue
        start_idx = module_idx_to_state_idx[module_idx]
        iteration_idx = next(state_idx
                             for state_idx in range(start_idx, len(states))
                             if states[state_idx].state_type == DetailedHMMStateType.CHOOSE_IF_ITERATE)
        if module.iterative_module:
            adj_list[iteration_idx][start_idx] = _make_edge(edge_type=DetailedHMMEdgeType.ITERATE_MODULE,
                                                                    log_prob=None,
                                                                    fst_module_idx=None)
        if module.iterative_gene and gene_intervals[module.gene_id][1] == module_idx:  # if the last module of the gene
            fst_module_in_gene_idx = gene_intervals[module.gene_id][0]
            edge_end = module_idx_to_state_idx[fst_module_in_gene_idx]
            adj_list[iteration_idx][edge_end] = _make_edge(edge_type=DetailedHMMEdgeType.ITERATE_GENE,
                                                           log_prob=None,
                                                           fst_module_idx=None)


def add_skip_genes(adj_list: Dict[int, Dict[int, DetailedHMMEdge]],
                   bgc_variant: BGC_Variant,
                   gene_intervals: Dict[GeneId, Tuple[int, int]],
                   module_idx_to_start_state_idx: Dict[int, int],
                   final_state_idx: int) -> None:
    _make_edge = partial(make_edge, bgc_variant=bgc_variant)

    # SKIP GENES AND FRAGMENTS IN THE MIDDLE
    for i, module in enumerate(bgc_variant.modules):
        module_start_state_idx = module_idx_to_start_state_idx[i]
        next_gene_start_state_idx = module_idx_to_start_state_idx.get(gene_intervals[module.gene_id][1] + 1, final_state_idx)
        if module.gene_id != bgc_variant.modules[i - 1].gene_id:  # module is the first in the gene
            adj_list[module_start_state_idx][next_gene_start_state_idx] = _make_edge(edge_type=DetailedHMMEdgeType.SKIP_GENE,
                                                                                     log_prob=None,
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
    _get_insert_at_start_emissions = partial(hmm_helper.get_insert_at_start_emissions,
                                             pks_domains_in_bgc=bgc_variant.has_pks_domains())

    # 1. Add inserts before modules which are at the start of fragments
    module_idx_to_insert_state_idx = {}
    for fragment_fst, _ in fragment_intervals.values():
        module = bgc_variant.modules[fragment_fst]
        states.append(DetailedHMMState(state_type=DetailedHMMStateType.INSERT_AT_START,
                                       emissions=_get_insert_at_start_emissions(module)))
        module_idx_to_insert_state_idx[fragment_fst] = len(states) - 1
        adj_list[len(states) - 1][module_idx_to_start_state_idx[fragment_fst]] = _make_edge(edge_type=DetailedHMMEdgeType.START_MATCHING,
                                                                                            log_prob=None,  # will be filled in later
                                                                                            fst_module_idx=None)
        adj_list[len(states) - 1][len(states) - 1] = _make_edge(edge_type=DetailedHMMEdgeType.INSERT,
                                                                log_prob=None,  # will be filled in later
                                                                fst_module_idx=None)

    # 2. No skips at the beginning
    adj_list[initial_state][module_idx_to_insert_state_idx[0]] = _make_edge(edge_type=DetailedHMMEdgeType.START_INSERTING_AT_START,
                                                                            log_prob=None,  # will be filled in later
                                                                            fst_module_idx=0)
    adj_list[initial_state][module_idx_to_start_state_idx[0]] = _make_edge(edge_type=DetailedHMMEdgeType.START_MATCHING,
                                                                           log_prob=None,
                                                                           fst_module_idx=0)

    # 3. Skip modules at the beginning
    module_idx_to_module_skip_state_idx = {}
    last_fragment_fst, _ = fragment_intervals[len(fragment_intervals) - 1]
    for i, module in enumerate(bgc_variant.modules[:last_fragment_fst]):  # skipping modules[:i+1]. Modules in the last fragment can't be skipped
        # if path comes to this state, the i-th module is skipped
        states.append(DetailedHMMState(state_type=DetailedHMMStateType.SKIP_MODULE_AT_START,
                                       emissions={}))
        module_idx_to_module_skip_state_idx[i] = len(states) - 1

        if bgc_variant.modules[i + 1].fragment_idx != bgc_variant.modules[i].fragment_idx:  # next module is the first in the fragment
            adj_list[len(states) - 1][module_idx_to_start_state_idx[i + 1]] = _make_edge(edge_type=DetailedHMMEdgeType.START_MATCHING,
                                                                                         log_prob=None,
                                                                                         fst_module_idx=None)
            adj_list[len(states) - 1][module_idx_to_insert_state_idx[i + 1]] = _make_edge(edge_type=DetailedHMMEdgeType.START_INSERTING_AT_START,
                                                                                          log_prob=None,
                                                                                          fst_module_idx=i + 1)
        if i > 0:
            prev_module_skip_state_idx = module_idx_to_module_skip_state_idx[i - 1]
            adj_list[prev_module_skip_state_idx][len(states) - 1] = _make_edge(edge_type=DetailedHMMEdgeType.SKIP_MODULE_AT_START,
                                                                               log_prob=None,
                                                                               fst_module_idx=i)  # skip the previous module

    if len(fragment_intervals) > 1:  # at least two fragments
        states.append(DetailedHMMState(state_type=DetailedHMMStateType.START_SKIP_MODULES_AT_START,
                                       emissions={}))
        adj_list[len(states) - 1][module_idx_to_module_skip_state_idx[0]] = _make_edge(edge_type=DetailedHMMEdgeType.SKIP_MODULE_AT_START,
                                                                                       log_prob=None,
                                                                                       fst_module_idx=0)
        adj_list[initial_state][len(states) - 1] = _make_edge(edge_type=DetailedHMMEdgeType.START_SKIP_MODULES_AT_START,
                                                              log_prob=None,
                                                              fst_module_idx=None)


def add_skip_at_the_end(states: List[DetailedHMMState],
                        adj_list: Dict[int, Dict[int, DetailedHMMEdge]],
                        bgc_variant: BGC_Variant,
                        module_idx_to_start_state_idx: Dict[int, int],
                        fragment_intervals: Dict[int, Tuple[int, int]],
                        hmm_helper: HMMHelper,
                        final_state_idx = -1) -> None:
    if len(fragment_intervals) == 1:  # bgc_variant has only one fragment
        return

    ET = DetailedHMMEdgeType
    ST = DetailedHMMStateType
    _make_edge = partial(make_edge, bgc_variant=bgc_variant)
    _get_insert_emissions = partial(hmm_helper.get_insert_emissions,
                                    pks_domains_in_bgc=bgc_variant.has_pks_domains())

    prev_skip_state_idx = None
    # SKIP AT THE END
    for module_idx, module in enumerate(bgc_variant.modules[fragment_intervals[1][0]:],
                                        start=fragment_intervals[1][0]):  # start from beginning of snd fragment because all fragments can't be skipped
        if module_idx < len(bgc_variant.modules) - 1:  # not the last module
            states.append(DetailedHMMState(state_type=ST.SKIP_MODULE_AT_END,
                                           emissions={}))
            skip_state_idx = len(states) - 1
        else:
            skip_state_idx = final_state_idx

        if prev_skip_state_idx is not None:  # connect to previous skip module state
            adj_list[prev_skip_state_idx][skip_state_idx] = _make_edge(edge_type=ET.SKIP_MODULE_AT_END,
                                                                       log_prob=None,
                                                                       fst_module_idx=module_idx)

        # start skipping if current module is the first in the fragment
        if bgc_variant.modules[module_idx - 1].fragment_idx != bgc_variant.modules[module_idx].fragment_idx:
            # note that emissions depend on the previous module: we're inserting *after it*
            states.extend([DetailedHMMState(state_type=ST.END_MATCHING,
                                            emissions={}),
                           DetailedHMMState(state_type=ST.INSERT_AT_END,
                                            emissions=_get_insert_emissions(bgc_variant.modules[module_idx - 1])),
                           DetailedHMMState(state_type=ST.END_INSERTING_AT_END,
                                            emissions={})])

            end_matching_idx, insert_at_end_idx, end_inserting_at_end = range(len(states) - 3, len(states))
            prev_module_start_idx = module_idx_to_start_state_idx[module_idx - 1]
            prev_module_match_idx = next(idx for idx in adj_list[prev_module_start_idx]
                                         if adj_list[prev_module_start_idx][idx].edge_type == ET.MATCH)

            adj_list[prev_module_match_idx][end_matching_idx] = _make_edge(edge_type=ET.END_MATCHING,
                                                                          log_prob=None,
                                                                          fst_module_idx=module_idx - 1)
            if module_idx > 1:  # just filter out an annoying case of a single-module first fragment: it doesn't break anything but skipping the whole BGC is pointless
                adj_list[prev_module_start_idx][end_matching_idx] = _make_edge(edge_type=ET.SKIP_MODULE_END_MATCHING,
                                                                               log_prob=None,
                                                                               fst_module_idx=module_idx - 1)

            adj_list[end_matching_idx][end_inserting_at_end] = _make_edge(edge_type=ET.NO_INSERTIONS,
                                                                          log_prob=None,
                                                                          fst_module_idx=None)
            adj_list[end_matching_idx][insert_at_end_idx] = _make_edge(edge_type=ET.START_INSERTING_AT_END,
                                                                       log_prob=None,
                                                                       fst_module_idx=module_idx - 1)
            adj_list[insert_at_end_idx][insert_at_end_idx] = _make_edge(edge_type=ET.INSERT_AT_END,
                                                                        log_prob=None,
                                                                        fst_module_idx=None)
            adj_list[insert_at_end_idx][end_inserting_at_end] = _make_edge(edge_type=ET.END_INSERTING,
                                                                           log_prob=None,
                                                                           fst_module_idx=None)
            adj_list[end_inserting_at_end][skip_state_idx] = _make_edge(edge_type=ET.SKIP_MODULE_AT_END,
                                                                           log_prob=None,
                                                                           fst_module_idx=module_idx)
        prev_skip_state_idx = skip_state_idx


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
    final_state_idx = -1  # at the very end will be replaced with len(states) - 1

    # ADD MODULES' STATES: MODULE_START -> MATCH -> INSERT_MONOMER -> MODULE_START -> ...
    module_idx_to_state_idx = add_modules(states,
                                          adj_list,
                                          bgc_variant,
                                          hmm_helper,
                                          gene_intervals,
                                          start_state_idx,
                                          final_state_idx)  # module_idx -> module_start_state_idx
    state_idx_to_module_idx = {v: k for k, v in module_idx_to_state_idx.items()}
    adj_list[start_state_idx][module_idx_to_state_idx[0]] = _make_edge(edge_type=DetailedHMMEdgeType.START_MATCHING,
                                                                       log_prob=None,  # will be set in the end
                                                                       fst_module_idx=0)


    # SKIP AT THE BEGINNING: SKIP_MODULE_AT_START -> INSERT_AT_STARTE -> // MODULE_START -> ...
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
                        fragment_intervals,
                        hmm_helper)

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
    # do it after topological sorting to enforce the order of modules
    '''
    add_skip_genes(adj_list,
                   bgc_variant,
                   gene_intervals,
                   module_idx_to_state_idx,
                   final_state_idx)
    '''

    # must be done after topological sorting because it creates cycles
    add_modules_and_genes_iterations(states,
                                     adj_list,
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
