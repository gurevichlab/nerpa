from __future__ import annotations
from typing import TYPE_CHECKING, Dict, Tuple
if TYPE_CHECKING:
    from src.hmm.hmm import HMM

from typing import List
from typing import Optional
from src.hmm.hmm_auxiliary_types import DetailedHMMStateType
from itertools import pairwise
import pulp
import shutil
from multiprocessing import Process
from pathlib import Path
from graphviz import Digraph


def render_with_fallback(graph: Digraph, out_png: Path, timeout: float = 10.0,
                         report_errors: bool = False) -> None:
    def clone_with_splines(s: str) -> Digraph:
        g = Digraph(engine=graph.engine, format=graph.format)
        g.body = list(graph.body)
        g.graph_attr.update(graph.graph_attr)
        g.node_attr.update(graph.node_attr)
        g.edge_attr.update(graph.edge_attr)
        g.attr(splines=s)  # <- important: goes through graphviz's attr machinery
        return g

    def render(s: str):
        clone_with_splines(s).render(out_png.with_suffix(""), cleanup=True)

    render("polyline")
    # p = Process(target=render, args=("curved",))
    # p.start()
    # p.join(timeout)

    # if p.is_alive():
    #     p.terminate()
    #     p.join()
    #     if report_errors:
    #         print(f"Warning: Graph rendering with 'curved' splines timed out after {timeout} seconds. "
    #               f"Falling back to 'polyline' splines, which may look worse.")
    #     render("polyline")


def create_nodes_layout(
        layers: List[List[int]],
        edges: List[Tuple[int, int]],
        labels: Dict[int, str],
        char_width: float = 0.6,
        label_padding: float = 0.8,
        edge_min_dx: float = 0.5,
        layer_gap: float = 0.6,
        layer_y_gap: float = 5,
) -> Dict[int, Tuple[float, float]]:
    # Nodes we care about: union of those in layers and those with labels
    nodes = sorted(node_idx for layer in layers for node_idx in layer)
    assert set(nodes) == set(labels.keys())

    def node_size(i: int) -> float:
        return len(labels[i]) * char_width + label_padding

    # layer index for y
    layer_of: Dict[int, int] = {}
    for li, layer in enumerate(layers):
        for n in layer:
            layer_of[n] = li

    # LP
    prob = pulp.LpProblem("layout", pulp.LpMinimize)
    x = {i: pulp.LpVariable(f"x_{i}", lowBound=None, upBound=None) for i in nodes}

    # Fix translation: x_0 = 0 (or smallest-index node)
    prob += (x[nodes[0]] == 0.0)

    # Objective: minimize x_max - x_min (here: last - first in index order)
    prob += x[nodes[-1]] - x[nodes[0]]

    # Index order constraints: x_{i+1} >= x_i + index_gap
    for u, v in edges:
        prob += x[v] >= x[u] + edge_min_dx

    # Same-layer non-overlap for adjacent-by-index nodes in that layer
    for layer in layers:
        for u, v in pairwise(layer):
            sep = (node_size(u) / 2.0) + layer_gap + (node_size(v) / 2.0)
            prob += x[v] >= x[u] + sep

    # Solve
    # In conda envs, CBC is usually provided as the `cbc` executable (coin-or-cbc).
    # PuLP's PULP_CBC_CMD expects a bundled CBC and may report "Not Available".
    cbc_path: str | None = shutil.which('cbc')
    status = prob.solve(pulp.COIN_CMD(path=cbc_path, msg=False))
    # status = prob.solve(pulp.PULP_CBC_CMD(msg=False))
    if pulp.LpStatus[status] != "Optimal":
        raise RuntimeError(f"LP did not solve optimally: {pulp.LpStatus[status]}")

    # Extract
    pos: Dict[int, Tuple[float, float]] = {}
    for i in nodes:
        xi = float(pulp.value(x[i]))
        yi = (len(layers) - layer_of[i]) * layer_y_gap
        pos[i] = (xi, yi)
    return pos

def state_idx_to_label(idx: int,
                       hmm: HMM,) -> str:
    ST = DetailedHMMStateType
    state_type = hmm.state_types[idx]
    if state_type == ST.MODULE_SUBGRAPH_ROOT:
        module_idx = hmm.module_start_states.index(idx)
        label = f'{idx}:{hmm.module_names[module_idx]}'
    else:
        label = f'{idx}:{state_type.name}'

    if len(label) > 10:
        label = label[:len(label) // 2] + '\n' + label[len(label) // 2:]

    return label


def draw_hmm(hmm: HMM,
             output_path: Optional[Path] = None,
             highlight_path: Optional[List[int]] = None,
             edge_weights: bool = True) -> Digraph:
    ST = DetailedHMMStateType

    graph = Digraph(format="svg", engine="neato")
    graph.graph_attr["neato"] = "-n2"
    graph.attr(
        overlap="false",
        splines="polyline",  # <- bends (most reliable)
        esep="+1",  # <- keep edges away from nodes (increase if needed)
        sep="+10",  # <- keep nodes separated (increase if needed)
    )
    #graph.attr(overlap="false", splines="true")  # nice routing; doesn’t move nodes

    # Define the layer types and reverse for top-to-bottom ordering
    layer_types = [[ST.SKIPPING_MODULES_AT_END],
                   [ST.INSERT_AT_END],
                   [ST.MATCHING_FINISHED,],
                   [ST.INSERT, ST.CHOOSE_IF_FINISH_MATCHING],
                   [ST.MATCH, ST.MATCH_POSSIBLE_ASSEMBLY_FINISH, ST.MATCH_LAST_MODULE, ST.MATCHING_CONTINUED],
                   [ST.INITIAL, ST.MODULE_SUBGRAPH_ROOT, ST.FINAL, ST.CHOOSE_IF_ITERATE_MODULE, ST.CHOOSE_IF_ITERATE_GENE, ST.CHOOSE_IF_ITERATE_BGC],
                   [ST.INSERT_AT_START],
                   [ST.CHOOSE_IF_START_MATCHING, ST.SKIPPING_MODULES_AT_START_FINISHED]]

    layers = [[] for _ in range(len(layer_types))]


    def state_idx_to_color(idx: int) -> str:
        if hmm.is_emitting_state(idx):
            if highlight_path is not None and idx in highlight_path:
                return "plum"
            else:
                return "lavenderblush"
        match hmm.state_types[idx]:
            case ST.MODULE_SUBGRAPH_ROOT:
                return "lightblue"
            case ST.INITIAL:
                return "grey"
            case ST.FINAL:
                return "grey"
            case _:
                return "white"

    # Assign nodes to their respective layers
    for state_idx in range(len(hmm.transitions)):
        for layer_idx, layer_state_types in enumerate(layer_types):
            if hmm.state_types[state_idx] in layer_state_types:
                layers[layer_idx].append(state_idx)

    labels = {i: state_idx_to_label(i, hmm)
              for i in range(len(hmm.transitions))}
    forward_edges = [
        (from_idx, to_idx)
        for from_idx, edges in enumerate(hmm.transitions)
        for to_idx, _lp in edges
        if to_idx > from_idx
    ]
    pos = create_nodes_layout(layers=layers,
                              labels=labels,
                              edges=forward_edges)  # (x,y) floats

    for state_idx in range(len(hmm.transitions)):
        x, y = pos[state_idx]
        graph.node(
            str(state_idx),
            label=labels[state_idx],
            shape="ellipse",
            fillcolor=state_idx_to_color(state_idx),
            style="filled",
            pos=f"{x},{y}!",  # the "!" pins it (critical)
        )

    path_edges = list(pairwise(highlight_path)) \
        if highlight_path is not None else []
    # Add actual edges between nodes
    for from_idx, edges in enumerate(hmm.transitions):
        for to_idx, lp in edges:
            edge_args = {
                'tail_name': str(from_idx),
                'head_name': str(to_idx),
                'arrowhead': "vee",
                'arrowsize': "1.5" if (from_idx, to_idx) in path_edges else "1",
                'color': "red" if (from_idx, to_idx) in path_edges else "black",
                'penwidth': "2" if (from_idx, to_idx) in path_edges else "1",
            }
            if edge_weights:
                edge_args['label'] = (
                    "-inf" if lp < -1e10 else
                    f"{lp:.2f}"
                )
            if from_idx == to_idx:
                match hmm.state_types[to_idx]:
                    case ST.INSERT:
                        edge_args['headport'] = "n"  # point upwards
                        edge_args['tailport'] = "n"
                    case ST.INSERT_AT_START | ST.INSERT_AT_END:
                        edge_args['headport'] = "w"  # point left
                        edge_args['tailport'] = "w"
                    case _:
                        raise ValueError(f"Unexpected self-loop at state {to_idx} of type {hmm.state_types[to_idx]}")

            graph.edge(**edge_args)

    # Optionally save the graph
    if output_path:
        render_with_fallback(graph, output_path, report_errors=True)
    return graph
