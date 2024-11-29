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
        pass  # defined after DetailedHMM definition

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
        """
        Draws the Hidden Markov Model (HMM) as a graph using Graphviz.

        Args:
            output_path (Optional[str]): The path to save the rendered graph.
                                         If None, the graph is returned without saving.
        Returns:
            Digraph: The Graphviz Digraph object representing the HMM.
        """
        graph = Digraph(format='png')
        graph.attr(rankdir='LR')  # Left-to-right graph orientation

        # Add nodes for all states without emissions
        for idx, state in enumerate(self.states):
            state_label = f"{idx}: {state.state_type.name}"
            graph.node(str(idx), label=state_label, shape="ellipse")

        # Add edges from adj_list
        for from_idx, edges in self.adj_list.items():
            for edge in edges:
                edge_label = edge.edge_type.name  # Only include edge type, no probabilities
                graph.edge(str(from_idx), str(edge.to), label=edge_label)

        # Optionally save the graph
        if output_path:
            graph.render(output_path, cleanup=True)
        return graph


DetailedHMM.from_bgc_variant = bgc_variant_to_detailed_hmm