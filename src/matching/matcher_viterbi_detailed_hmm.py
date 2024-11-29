from __future__ import annotations
from typing import List, Tuple, NamedTuple, Dict, Optional
from dataclasses import dataclass
from enum import Enum, auto

from src.data_types import NRP_Monomer, BGC_Module, BGC_Variant, GeneId, LogProb
from src.matching.matcher_viterbi_types import DetailedHMMEdgeType, DetailedHMMStateType, DetailedHMMState, DetailedHMMEdge, HMM
from collections import defaultdict
from src.matching.bgc_to_hmm import bgc_variant_to_detailed_hmm
from graphviz import Digraph
from pathlib import Path


@dataclass
class DetailedHMM:
    states: List[DetailedHMMState]
    adj_list: Dict[int, List[DetailedHMMEdge]]
    bgc_variant: BGC_Variant
    state_idx_to_module_idx: Dict[int, int]

    @classmethod
    def from_bgc_variant(cls, bgc_variant: BGC_Variant) -> DetailedHMM:
        return bgc_variant_to_detailed_hmm(DetailedHMM, bgc_variant)

    def to_hmm(self) -> HMM:
        num_states = len(self.states)
        adj_list = [[(edge.to, edge.log_prob) for edge in self.adj_list[u]]
                    for u in range(num_states)]

        emission_log_probs = [[state.emissions[mon]
                               for mon in sorted(state.emissions, key=lambda m: int(m))]
                              for state in self.states]
        return HMM(adj_list=adj_list,
                   emission_log_probs=emission_log_probs)

    def draw(self, output_path: Optional[Path] = None) -> Digraph:
        graph = Digraph(format='png')
        graph.attr(rankdir='BT')  # Top-to-bottom layout

        # Separate nodes into layers
        bottom_layer = []
        middle_layer = []
        top_layer = []

        nodes_to_draw = [idx for idx, state in enumerate(self.states)
                         if state.state_type in {
                             DetailedHMMStateType.INITIAL,
                             DetailedHMMStateType.MODULE_START,
                             DetailedHMMStateType.FINAL,
                            DetailedHMMStateType.MATCH,
                            DetailedHMMStateType.INSERT_MONOMER}]
        for idx, state in enumerate(self.states):
            if state.state_type in {
                DetailedHMMStateType.INITIAL,
                DetailedHMMStateType.MODULE_START,
                DetailedHMMStateType.FINAL,
            }:
                if state.state_type == DetailedHMMStateType.MODULE_START:
                    module = self.bgc_variant.modules[self.state_idx_to_module_idx[idx]]
                    label = f'{idx}:{module.gene_id}:{module.a_domain_idx}'
                else:
                    label = f'{idx}:{state.state_type.name}'
                bottom_layer.append((str(idx), label))
            elif state.state_type == DetailedHMMStateType.MATCH:
                middle_layer.append((str(idx), f'{idx}:{state.state_type.name}'))
            elif state.state_type == DetailedHMMStateType.INSERT_MONOMER:
                top_layer.append((str(idx), f'{idx}:{state.state_type.name}'))

        # Define layers
        with graph.subgraph() as sub:
            sub.attr(rank="same")
            for node, label in bottom_layer:
                sub.node(node, label=label, shape="ellipse")

        with graph.subgraph() as sub:
            sub.attr(rank="same")
            for node, label in middle_layer:
                sub.node(node, label=label, shape="ellipse")

        with graph.subgraph() as sub:
            sub.attr(rank="same")
            for node, label in top_layer:
                sub.node(node, label=label, shape="ellipse")

        # Add edges
        for from_idx, edges in self.adj_list.items():
            for edge in edges:
                if from_idx not in nodes_to_draw or edge.to not in nodes_to_draw:
                    continue
                if str(from_idx) == str(edge.to):  # Self-loop case
                    graph.edge(
                        str(from_idx),
                        str(edge.to),
                        label='',#edge.edge_type.name,
                        arrowhead="vee",
                        arrowsize="0.5",
                        headport="n",  # Anchor at the top of the node
                        tailport="n",  # Start at the top of the node
                    )
                else:
                    graph.edge(
                        str(from_idx),
                        str(edge.to),
                        label='',#edge.edge_type.name,
                        arrowhead="vee",
                        arrowsize="0.5",
                    )

        # Optionally save the graph
        if output_path:
            if '.' in output_path.name and not output_path.name.endswith('.png'):
                graph.render(output_path.name.rsplit('.', 1)[0], cleanup=True)
            else:
                graph.render(output_path.name, cleanup=True)
        return graph

