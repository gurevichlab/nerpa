from typing import List, NamedTuple, Union
import networkx as nx
import copy
from itertools import chain, permutations


def hamiltonian_path(G: nx.DiGraph,
                     source: int) -> Union[List[int], None]:
    ''' finds a hamiltonian path in G starting at source in an exhaustive manner '''
    visited = set()

    def dfs(u: int) -> Union[List[int], None]:
        visited.add(u)
        if len(visited) == len(G.nodes()):
            return [u]

        for v in G.successors(u):
            if v not in visited and (path := dfs(v)):
                return path + [u]

        visited.remove(u)
        return None

    if path := dfs(source):
        return path[::-1]
    else:
        return None


def parse_as_simple_cycle(G: nx.DiGraph) -> Union[List[int], None]:
    try:
        cycle_edges = nx.find_cycle(G)
        if len(cycle_edges) == len(G.edges()):  # since all edges are different, this implies the sets are equal as well
            return [u for u, v in cycle_edges]
    except nx.NetworkXNoCycle:
        return None


class BackboneSequence(NamedTuple):
    node_idxs: List[int]
    is_cyclic: bool = False


def putative_backbones(G_: nx.DiGraph, min_nodes: int = 2) -> List[BackboneSequence]:
    ''' Tries to parse each component of G either as a simple cycle or as a hamiltonian path '''
    G = copy.deepcopy(G_)  # there's no need for this but modifying an argument makes me anxious
    breakage_points = [node for node in G.nodes()
                       if G.in_degree(node) == 0 and G.out_degree(node) > 1 \
                       or G.in_degree(node) > 1 and G.out_degree(node) == 0]  # sinks and sources of degree > 1
    G.remove_nodes_from(breakage_points)

    backbone_sequences = []
    for component in nx.connected_components(nx.Graph(G)):
        if len(component) < min_nodes:
            continue

        Gs = G.subgraph(component)
        if cycle := parse_as_simple_cycle(Gs):
            backbone_sequences.append(BackboneSequence(node_idxs=cycle,
                                                       is_cyclic=True))
        elif ham_path := next(filter(None, (hamiltonian_path(Gs, u) for u in Gs)), None):
            backbone_sequences.append(BackboneSequence(node_idxs=ham_path))

    return backbone_sequences


def permuted_backbones(backbones: List[BackboneSequence]) -> List[BackboneSequence]:
    if 2 <= len(backbones) <= 3 and not any(backbone.is_cyclic for backbone in backbones):
        joined_paths = (list(chain(*(backbone.node_idxs for backbone in perm_backbones)))
                        for perm_backbones in permutations(backbones))
        return [BackboneSequence(node_idxs=joined_idxs,
                                 is_cyclic=False)
                for joined_idxs in joined_paths]
    else:
        return []
