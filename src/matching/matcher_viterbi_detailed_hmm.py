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
                       [DetailedHMMStateType.INSERT],
                       [DetailedHMMStateType.MATCH],
                       [DetailedHMMStateType.INITIAL, DetailedHMMStateType.MODULE_START, DetailedHMMStateType.FINAL],
                       [DetailedHMMStateType.INSERT_AT_START],
                       [DetailedHMMStateType.SKIP_MODULE_AT_START],
                       [DetailedHMMStateType.SKIP_GENE_AT_START],
                       [DetailedHMMStateType.SKIP_FRAGMENT_AT_START]][::-1]  # Reverse to match top-to-bottom order

        layers = [[] for _ in range(len(layer_types))]

        def state_idx_to_label(idx: int) -> str:
            state = self.states[idx]
            if state.state_type == DetailedHMMStateType.MODULE_START:
                module = self.bgc_variant.modules[self.state_idx_to_module_idx[idx]]
                return f'{idx}:F{module.fragment_idx}:{module.gene_id}:{module.a_domain_idx}'
            else:
                return f'{idx}:{state.state_type.name}'

        def state_idx_to_color(idx: int) -> str:
            match self.states[idx].state_type:
                case DetailedHMMStateType.MODULE_START:
                    return "lightblue"
                case DetailedHMMStateType.INITIAL:
                    return "grey"
                case DetailedHMMStateType.FINAL:
                    return "grey"
                case _:
                    return "white"

        # Assign nodes to their respective layers
        for idx, state in enumerate(self.states):
            for layer_idx, layer_state_types in enumerate(layer_types):
                if state.state_type in layer_state_types:
                    layers[layer_idx].append(idx)

        # Sort nodes within each layer to enforce order
        for layer in layers:
            layer.sort()

        dummy_nodes = [f'dummy{i}' for i in range(len(layer_types))]  # Dummy nodes for enforcing horizontal order
        # Define layers and enforce horizontal order using invisible chains
        for layer_idx, (layer, dummy_node) in enumerate(zip(layers, dummy_nodes)):
            with graph.subgraph() as sub:
                sub.attr(rank="same")
                sub.node(dummy_node, label="", shape="none", width="0", height="0", style="invis")  # Add dummy node

                # Create an invisible chain for horizontal order
                invisible_chain = []
                for i in range(len(layer)):
                    chain_node = f'chain_{layer_idx}_{i}'
                    invisible_chain.append(chain_node)
                    sub.node(chain_node, label="", shape="none", width="0", height="0",
                             style="invis")  # Invisible chain node

                # Connect the invisible chain sequentially
                for chain1, chain2 in zip(invisible_chain[:-1], invisible_chain[1:]):
                    sub.edge(chain1, chain2, style="invis")

                # Anchor each actual node to its position in the invisible chain
                for state_idx, chain_node in zip(layer, invisible_chain):
                    sub.node(str(state_idx), label=state_idx_to_label(state_idx), shape="ellipse",
                             fillcolor=state_idx_to_color(state_idx),
                             style="filled")
                    sub.edge(chain_node, str(state_idx), style="invis")
                    sub.edge(dummy_node, str(state_idx), style="invis")  # Tie to dummy node

        # Connect dummy nodes **sequentially** to enforce correct layer order
        for dummy1, dummy2 in zip(dummy_nodes[:-1], dummy_nodes[1:]):
            graph.edge(dummy1, dummy2, style="invis")  # Connect only adjacent dummy nodes

        # Add actual edges between nodes
        for from_idx, edges in self.adj_list.items():
            for edge in edges:
                edge_args = {
                    'tail_name': str(from_idx),
                    'head_name': str(edge.to),
                    'label':  edge.edge_type.name,
                    'arrowhead': "vee",
                    'arrowsize': "0.5"
                }
                if str(from_idx) == str(edge.to):  # Self-loop case
                    if self.states[from_idx].state_type == DetailedHMMStateType.INSERT:
                        edge_args['headport'] = "n"  # point upwards
                        edge_args['tailport'] = "n"
                    else:  # INSERT_AT_START
                        edge_args['headport'] = "s"  # point downwards
                        edge_args['tailport'] = "s"
                graph.edge(**edge_args)

        # Optionally save the graph
        if output_path:
            if '.' in output_path.name and not output_path.name.endswith('.png'):
                graph.render(output_path.name.rsplit('.', 1)[0], cleanup=True)
            else:
                graph.render(output_path.name, cleanup=True)
        return graph
