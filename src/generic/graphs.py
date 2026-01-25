from typing import (
    Dict,
    Iterable,
    List,
    NamedTuple,
    Optional,
    TypeVar,
    Tuple,
    Union, Callable, Hashable, Generic
)
import networkx as nx
import copy
from itertools import chain, permutations, pairwise

from networkx.algorithms.isomorphism.isomorphvf2 import DiGraphMatcher


def hamiltonian_path(G: nx.DiGraph,
                     source: int) -> Optional[List[int]]:
    ''' finds a hamiltonian path in G starting at source in an exhaustive manner '''
    visited = set()

    def dfs(u: int) -> Optional[List[int]]:
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


def putative_backbones(G_: nx.DiGraph, min_nodes: int = None) -> List[BackboneSequence]:
    ''' Tries to parse each component of G either as a simple cycle or as a hamiltonian path '''
    G = copy.deepcopy(G_)  # there's no need for this but modifying an argument makes me anxious
    breakage_points = [node for node in G.nodes()
                       if G.in_degree(node) == 0 and G.out_degree(node) > 1 \
                       or G.in_degree(node) > 1 and G.out_degree(node) == 0]  # sinks and sources of degree > 1
    G.remove_nodes_from(breakage_points)

    if min_nodes is None:
        backbone_sequences = [BackboneSequence(node_idxs=[breakage_point], is_cyclic=False)
                              for breakage_point in breakage_points]
    else:
        backbone_sequences = []

    for component in nx.connected_components(nx.Graph(G)):
        if min_nodes is not None and len(component) < min_nodes:
            continue

        Gs = G.subgraph(component)
        if cycle := parse_as_simple_cycle(Gs):
            backbone_sequences.append(BackboneSequence(node_idxs=cycle,
                                                       is_cyclic=True))
        else:
            ham_paths = list(filter(None, (hamiltonian_path(Gs, u) for u in Gs)))
            if len(ham_paths) > 1:
                print('WARNING! Multiple hamiltonian paths found. Proceeding with the first one')
                #raise ValueError('Multiple hamiltonian paths found')  # TODO: handle this case
            if ham_paths:
                backbone_sequences.append(BackboneSequence(node_idxs=ham_paths[0]))

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


NodeId = TypeVar('NodeId')
Edge = Tuple[NodeId, NodeId]


def get_interior_edges(nodes: List[NodeId],
                       edges: Iterable[Edge]) -> List[Edge]:
    return [(u, v)
            for u, v in edges
            if u in nodes and v in nodes]


def get_boundary_edges(nodes: List[NodeId],
                       edges: Iterable[Edge]) -> List[Edge]:
    return [(u, v)
            for u, v in edges
            if any([u in nodes and v not in nodes,
                    v in nodes and u not in nodes])]


def shortest_path_through(G: nx.DiGraph,
                          include_nodes: Optional[List[NodeId]] = None,
                          exclude_nodes: Optional[List[NodeId]] = None,
                          edge_weights: Dict[Tuple[NodeId, NodeId], float] = None) -> Optional[List[NodeId]]:
    '''
    finds the shortest path that goes through include_nodes and avoids exclude_nodes
    the path goes through the include_nodes in the order they are given
    include_nodes can have repetitions, in which case the path will go through them multiple times
    '''

    if include_nodes is None:
        include_nodes = []
    if exclude_nodes is None:
        exclude_nodes = []
    if edge_weights is None:
        edge_weights = {edge: 1 for edge in G.edges()}
    _G = G.copy()
    _G.remove_nodes_from(exclude_nodes)

    edge_weight_fun = lambda u, v, d: edge_weights[(u, v)]
    path = [include_nodes[0]]
    for start, finish in pairwise(include_nodes):
        path_segment = nx.shortest_path(_G, start, finish, weight=edge_weight_fun)
        if path_segment is None:
            return None
        path.extend(path_segment[1:])
    return path

NodeLabel = TypeVar('NodeLabel')

def graphs_one_substitution_away(G1: nx.DiGraph,
                                 G2: nx.DiGraph,
                                 nodes_comparator: Callable[[NodeLabel, NodeLabel], bool],
                                 label_key: str = 'label') -> bool:
    ''' checks if G1 can be transformed into G2 by substituting a single node '''
    
    # Quick structural filters
    if G1.number_of_nodes() != G2.number_of_nodes():
        return False
    if G1.number_of_edges() != G2.number_of_edges():
        return False

    # Optional: cheap degree-sequence filter (often helps a lot)
    if sorted(dict(G1.in_degree()).values()) != sorted(dict(G2.in_degree()).values()):
        return False
    if sorted(dict(G1.out_degree()).values()) != sorted(dict(G2.out_degree()).values()):
        return False

    # Enumerate all directed graph isomorphisms, ignoring labels
    matcher = DiGraphMatcher(G1, G2)

    # If there are no isomorphisms structurally, we're done
    return any(
        sum(not nodes_comparator(G1.nodes[u1].get(label_key),
                                 G2.nodes[u2].get(label_key))
            for u1, u2 in mapping.items()) == 1
        for mapping in matcher.isomorphisms_iter()
    )

T = TypeVar("T", bound=Hashable)

class DSU(Generic[T]):
    def __init__(self) -> None:
        self.parent: Dict[T, T] = {}

    def find(self, x: T) -> T:
        if x not in self.parent:
            self.parent[x] = x
            return x

        p = self.parent[x]
        if p != x:
            self.parent[x] = self.find(p)
        return self.parent[x]

    def union(self, a: T, b: T) -> T:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return ra

        root = min(ra, rb)   # requires orderable T
        self.parent[ra] = root
        self.parent[rb] = root
        return root

    # Iterate over all elements ever seen
    def __iter__(self) -> Iterable[T]:
        return iter(self.parent)

    # Optional convenience
    def items(self) -> Dict[T, T]:
        return {x: self.find(x) for x in self.parent}