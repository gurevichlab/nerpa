from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Dict, Optional

from src.antismash_parsing.bgc_variant_types import BGC_Variant, GeneId, BGC_Module
from src.general_type_aliases import LogProb
from src.hmm.hmm_constructor.hmm_constructor_iteration import add_iteration_subgraph
from src.hmm.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    DetailedHMMStateType,
    DetailedHMMState,
    DetailedHMMEdge,
    StateIdx,
)
from src.hmm.hmm_constructor.genes_fragments_intervals import get_genes_intervals
from collections import defaultdict

from src.hmm.hmm_constructor.hmm_constructor_module_subgraphs import add_module_subgraphs
from src.hmm.hmm_constructor.hmm_constructor_skips_at_end import add_skips_at_end
from src.hmm.hmm_constructor.hmm_constructor_skips_at_start import add_skips_at_start
from src.hmm.hmm_constructor.hmm_constructor_state_edge_context_relations import STATE_TYPE_TO_EDGE_TYPES, MATCHING_STATE_TYPES
from src.hmm.hmm_constructor.hmm_contructor_make_edge import make_edge
from src.hmm.hmm_scoring_helper import HMMHelper
import networkx as nx

from src.monomer_names_helper import NRP_Monomer

ET = DetailedHMMEdgeType
ST = DetailedHMMStateType

@dataclass
class HMM_Construction_Config:
    ASSEMBLY_LINE_BORDERS_ONLY_AT_MARKERS: bool = True

class HMM_Constructor:
    bgc_variant: BGC_Variant
    hmm_helper: HMMHelper
    construction_cfg: HMM_Construction_Config

    gene_borders: Dict[GeneId, Tuple[int, int]]

    states: List[DetailedHMMState] = None
    adj_list: Dict[StateIdx, Dict[StateIdx, DetailedHMMEdge]] = None
    module_idx_to_subgraph_root: Dict[int, StateIdx] = None
    initial_state: StateIdx = None
    final_state: StateIdx = None

    def __init__(self,
                 bgc_variant: BGC_Variant,
                 hmm_helper: HMMHelper):
        self.bgc_variant = bgc_variant
        self.hmm_helper = hmm_helper
        self.gene_borders = get_genes_intervals(bgc_variant.modules)

        # TODO: put into config
        self.construction_cfg = HMM_Construction_Config(
            ASSEMBLY_LINE_BORDERS_ONLY_AT_MARKERS=False
        )

    def make_edge(self,
                  state: DetailedHMMState,
                  edge_type: DetailedHMMEdgeType) -> DetailedHMMEdge:
        return make_edge(helper=self, state=state, edge_type=edge_type)

    def add_edge(self,
                 from_idx: StateIdx,
                 to_idx: StateIdx,
                 edge_type: DetailedHMMEdgeType):
        self.adj_list[from_idx][to_idx] = self.make_edge(state=self.states[from_idx],
                                                         edge_type=edge_type)

    def get_emissions(self,
                      module: BGC_Module) -> Dict[NRP_Monomer, LogProb]:
        return self.hmm_helper.get_emissions(module,
                                             pks_domains_in_bgc=self.bgc_variant.has_pks_domains())

    def get_insert_emissions(self,
                             module: BGC_Module) -> Dict[NRP_Monomer, LogProb]:
        return self.hmm_helper.get_insert_emissions(module,
                                                    pks_domains_in_bgc=self.bgc_variant.has_pks_domains())

    def add_iteration_subgraph(self,
                               module_idx: int,
                               node_after_subgraph: StateIdx) -> StateIdx:
        return add_iteration_subgraph(constructor=self,
                                      module_idx=module_idx,
                                      node_after_subgraph=node_after_subgraph)

    def topsort_states(self):
        g = nx.DiGraph()
        for from_idx, edges in self.adj_list.items():
            for to_idx, edge in edges.items():
                if (from_idx != to_idx and
                        edge.edge_type not in (ET.ITERATE_MODULE, ET.ITERATE_GENE)):
                    g.add_edge(from_idx, to_idx)

        topsorted_state_indices = list(nx.topological_sort(g))
        self.states = [self.states[i] for i in topsorted_state_indices]
        old_to_new_idx = {old_idx: new_idx for new_idx, old_idx in enumerate(topsorted_state_indices)}
        new_adj_list = defaultdict(dict)
        for old_from_idx, edges in self.adj_list.items():
            new_from_idx = old_to_new_idx[old_from_idx]
            for old_to_idx, edge in edges.items():
                new_to_idx = old_to_new_idx[old_to_idx]
                new_adj_list[new_from_idx][new_to_idx] = edge

        self.adj_list = new_adj_list
        self.initial_state = old_to_new_idx[self.initial_state]
        self.final_state = old_to_new_idx[self.final_state]
        self.module_idx_to_subgraph_root = {
            module_idx: old_to_new_idx[old_state_idx]
            for module_idx, old_state_idx in self.module_idx_to_subgraph_root.items()
        }

    def add_final_state(self):
        old_final_state = self.final_state
        self.states.append(DetailedHMMState(state_type=ST.FINAL))
        self.final_state = len(self.states) - 1
        for from_idx, edges_dict in self.adj_list.items():
            if old_final_state in edges_dict:
                edge = edges_dict[old_final_state]
                self.adj_list[from_idx][self.final_state] = edge
                del self.adj_list[from_idx][old_final_state]

    def check_hmm_sanity(self):
        assert self.initial_state == 0, "Initial state index should be 0"
        assert self.final_state == len(self.states) - 1, "Final state should be the last state"
        assert sorted(self.adj_list.keys()) == list(range(len(self.states) - 1)), "Some states are missing in the adjacency list"  # -1 because final state has no outgoing edges
        assert all(to_idx >= from_idx
                   for from_idx, edges_dict in self.adj_list.items()
                   for to_idx, edge in edges_dict.items()
                     if edge.edge_type not in (ET.ITERATE_MODULE,
                                                ET.ITERATE_GENE)), "HMM states are not properly topsorted"
        assert all(len(edges_dict) <= 2
                   for edges_dict in self.adj_list.values()), "There should be at most 2 outgoing edges per state"
        for state_idx, edges_dict in self.adj_list.items():
            edges_types = {edge.edge_type for edge in edges_dict.values()}
            expected_edges_types = STATE_TYPE_TO_EDGE_TYPES[self.states[state_idx].state_type]
            assert edges_types == expected_edges_types, \
                (f"Wrong edge types for state type {self.states[state_idx].state_type}\n"
                 f"Expected: {expected_edges_types}"
                 f"Got: {edges_types}")

    def convert_to_detailed_hmm(self, cls) -> 'DetailedHMM':
        state_idx_to_module_idx: Dict[StateIdx, Optional[int]] = {
            state_idx: state.related_module_idx
            for state_idx, state in enumerate(self.states)
        }
        match_states: Dict[int, StateIdx] = {
            state.related_module_idx: state_idx
            for state_idx, state in enumerate(self.states)
            if state.state_type in MATCHING_STATE_TYPES
        }

        match_states_list = [match_states[module_idx] for module_idx in range(len(self.bgc_variant.modules))]
        module_idx_to_subgraph_root_list = [self.module_idx_to_subgraph_root[module_idx]
                                           for module_idx in range(len(self.bgc_variant.modules))]

        return cls(states=self.states,
                   transitions=self.adj_list,
                   bgc_variant=self.bgc_variant,
                   state_idx_to_module_idx=state_idx_to_module_idx,
                   _module_idx_to_state_idx=module_idx_to_subgraph_root_list,
                   _module_idx_to_match_state_idx=match_states_list,
                   start_state_idx=self.initial_state,
                   final_state_idx=self.final_state,
                   hmm_helper=self.hmm_helper)

    def build_hmm(self, cls, debug=False):  # cls = DetailedHMM, have to avoid circular import
        debug = False
        self.states = [DetailedHMMState(state_type=ST.INITIAL,
                                        related_module_idx=0)]
        self.initial_state = 0
        self.final_state = -1  # will be updated at step 2
        self.adj_list = defaultdict(dict)

        # 0. Plant subgraph roots
        self.module_idx_to_subgraph_root = {i: i + 1  # +1 to reserve 0 for the initial state
                                            for i in range(len(self.bgc_variant.modules))}
        self.states.extend(
            DetailedHMMState(state_type=ST.MODULE_SUBGRAPH_ROOT,
                             related_module_idx=i)
            for i in range(len(self.bgc_variant.modules))
        )

        # 1. Build preliminary HMM
        add_module_subgraphs(self)
        if debug:
            self.convert_to_detailed_hmm(cls).draw(Path("added_module_subgraphs.png"))  # debugging
        add_skips_at_start(self, cls)
        if debug:
            self.convert_to_detailed_hmm(cls).draw(Path("added_skips_at_start.png"))
        add_skips_at_end(self)
        if debug:
            self.convert_to_detailed_hmm(cls).draw(Path("added_skips_at_end.png"))

        # 2. Add final state and redirect edges from the stub to it
        self.add_final_state()

        # 3. Topsort states (for Viterbi algorithm)
        self.topsort_states()
        if debug:
            self.convert_to_detailed_hmm(cls).draw(Path("hmm_built.png"))

        # 4. Sanity check
        self.check_hmm_sanity()

        return self.convert_to_detailed_hmm(cls)



