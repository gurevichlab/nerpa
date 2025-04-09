from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from src.matching.detailed_hmm import DetailedHMM
from typing import List
from pathlib import Path
from typing import Optional
from graphviz import Digraph
from src.matching.hmm_auxiliary_types import DetailedHMMStateType, DetailedHMMState, DetailedHMMEdge, \
    DetailedHMMEdgeType
from itertools import pairwise


def draw_hmm(hmm: DetailedHMM,
             output_path: Optional[Path] = None,
             highlight_path: Optional[List[int]] = None,
             edge_labels: bool = False) -> Digraph:
    ST = DetailedHMMStateType
    ET = DetailedHMMEdgeType

    graph = Digraph(format='png')
    graph.attr(rankdir='BT')  # Bottom-to-top layout to match desired top-to-bottom order

    # Define the layer types and reverse for top-to-bottom ordering
    layer_types = [[ST.SKIP_MODULE_AT_END],
                   [ST.END_INSERTING_AT_END],
                   [ST.INSERT_AT_END],
                   [ST.INSERT, ST.END_MATCHING],
                   [ST.MATCH],
                   [ST.INITIAL, ST.MODULE_START, ST.FINAL],
                   [ST.INSERT_AT_START],
                   [ST.SKIP_MODULE_AT_START]][::-1]  # Reverse to match top-to-bottom order

    layers = [[] for _ in range(len(layer_types))]

    def state_idx_to_label(idx: int) -> str:
        state = hmm.states[idx]
        if state.state_type == ST.MODULE_START:
            module = hmm.bgc_variant.modules[hmm.state_idx_to_module_idx[idx]]
            return f'{idx}:F{module.fragment_idx}:{module.gene_id}:{module.a_domain_idx}'
        else:
            return f'{idx}:{state.state_type.name}'

    def state_idx_to_color(idx: int) -> str:
        match hmm.states[idx].state_type:
            case ST.MODULE_START:
                return "lightblue"
            case ST.INITIAL:
                return "grey"
            case ST.FINAL:
                return "grey"
            case _:
                return "white"

    # Assign nodes to their respective layers
    for idx, state in enumerate(hmm.states):
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

    path_edges = list(pairwise(highlight_path)) \
        if highlight_path is not None else []
    # Add actual edges between nodes
    for from_idx, edges_dict in hmm.transitions.items():
        for to_idx, edge in edges_dict.items():
            edge_args = {
                'tail_name': str(from_idx),
                'head_name': str(to_idx),
                'arrowhead': "vee",
                'arrowsize': "1.5" if (from_idx, to_idx) in path_edges else "1",
                'color': "red" if (from_idx, to_idx) in path_edges else "black",
                'penwidth': "2" if (from_idx, to_idx) in path_edges else "1"
            }
            if edge_labels:
                edge_args['label'] = edge.edge_type.name
            if edge.edge_type == ET.INSERT:  # Self-loop case
                edge_args['headport'] = "n"  # point upwards
                edge_args['tailport'] = "n"
            if edge.edge_type == ET.INSERT_AT_START:
                edge_args['headport'] = "s"  # point downwards
                edge_args['tailport'] = "s"

            if (hmm.states[from_idx].state_type == ST.SKIP_MODULE_AT_END and
                    edge.edge_type == ET.SKIP_MODULE_AT_END):
                edge_args['headport'] = "n"  # point upwards
                edge_args['tailport'] = "n"
            graph.edge(**edge_args)

    # Optionally save the graph
    if output_path:
        if output_path.name.endswith('.png'):
            fname = str(output_path).rsplit('.', 1)[0]
            graph.render(fname, cleanup=True)
    return graph
