from __future__ import annotations
from typing import List, Tuple, NamedTuple, Dict, Optional
from dataclasses import dataclass
from enum import Enum, auto

from src.data_types import NRP_Monomer, BGC_Module, BGC_Variant, GeneId, LogProb
from src.matching.matcher_viterbi_types import DetailedHMMEdgeType, DetailedHMMStateType, DetailedHMMState, DetailedHMMEdge, HMM
from collections import defaultdict
from src.matching.bgc_to_hmm import bgc_variant_to_detailed_hmm
from itertools import pairwise
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
        graph.attr(rankdir='BT')  # Bottom-to-top layout to match desired top-to-bottom order

        # Define the layer types and reverse for top-to-bottom ordering
        layer_types = [[DetailedHMMStateType.SKIP_FRAGMENT_END],
                       [DetailedHMMStateType.INSERT_MONOMER],
                       [DetailedHMMStateType.MATCH],
                       [DetailedHMMStateType.INITIAL, DetailedHMMStateType.MODULE_START, DetailedHMMStateType.FINAL],
                       [DetailedHMMStateType.INSERT_AT_START],
                       [DetailedHMMStateType.SKIP_MODULE_AT_START],
                       [DetailedHMMStateType.SKIP_GENE_AT_START],
                       [DetailedHMMStateType.SKIP_FRAGMENT_AT_START]][::-1]  # Reverse to match top-to-bottom order

        layers = [[] for _ in range(len(layer_types))]
        dummy_nodes = [f'dummy{i}' for i in range(len(layer_types))]

        # Assign nodes to their respective layers
        for idx, state in enumerate(self.states):
            if state.state_type == DetailedHMMStateType.MODULE_START:
                module = self.bgc_variant.modules[self.state_idx_to_module_idx[idx]]
                label = f'{idx}:F{module.fragment_idx}:{module.gene_id}:{module.a_domain_idx}'
            else:
                label = f'{idx}:{state.state_type.name}'

            for layer_idx, layer_state_types in enumerate(layer_types):
                if state.state_type in layer_state_types:
                    layers[layer_idx].append((str(idx), label))
                    break

        # Define layers and connect dummy nodes to enforce rank
        for layer_idx, (layer, dummy_node) in enumerate(zip(layers, dummy_nodes)):
            with graph.subgraph() as sub:
                sub.attr(rank="same")
                sub.node(dummy_node, label="", shape="none", width="0", height="0", style="invis")  # Add dummy node

                if layer_idx == 4:  # Make Layer 4 more sparse
                    for i, (node, label) in enumerate(layer):
                        # Add the actual node
                        sub.node(node, label=label, shape="ellipse",
                                 fillcolor="lightblue" if sum(c == ':' for c in label) > 1 else "grey",
                                 style="filled")
                        if i > 0:
                            dummy_spacing = f'spacer_{layer_idx}_{i-1}'  # previous dummy
                            sub.edge(dummy_spacing, node, style="invis")
                        # Add an invisible dummy node for spacing
                        if i < len(layer) - 1:  # No dummy node after the last node
                            dummy_spacing = f'spacer_{layer_idx}_{i}'
                            sub.node(dummy_spacing, label="", shape="none", width="0", height="0", style="invis")
                            sub.edge(node, dummy_spacing, style="invis")
                else:
                    # Regular layer handling for all other layers
                    for node, label in layer:
                        sub.node(node, label=label, shape="ellipse")
                    for node, label in layer:
                        graph.edge(dummy_node, node, style="invis")  # Attach all layer nodes to the dummy node

        # Connect dummy nodes **sequentially** to enforce correct layer order
        for dummy1, dummy2 in zip(dummy_nodes[:-1], dummy_nodes[1:]):
            graph.edge(dummy1, dummy2, style="invis")  # Connect only adjacent dummy nodes

        # Add actual edges between nodes
        for from_idx, edges in self.adj_list.items():
            for edge in edges:
                if str(from_idx) == str(edge.to):  # Self-loop case
                    graph.edge(
                        str(from_idx),
                        str(edge.to),
                        label='',  # Optional: edge.edge_type.name
                        arrowhead="vee",
                        arrowsize="0.5",
                        headport="n",  # Anchor at the top of the node
                        tailport="n",  # Start at the top of the node
                    )
                else:
                    graph.edge(
                        str(from_idx),
                        str(edge.to),
                        label='',  # Optional: edge.edge_type.name
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
