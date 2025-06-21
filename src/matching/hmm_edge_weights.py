from itertools import pairwise, islice
from pathlib import Path
from typing import Dict, Tuple, TYPE_CHECKING, Union

from src.data_types import Prob, LogProb
from src.matching.hmm_scoring_config import HMMScoringConfig

if TYPE_CHECKING:
    from src.matching.detailed_hmm import DetailedHMM

from src.matching.hmm_auxiliary_types import DetailedHMMEdgeType, DetailedHMMEdge, DetailedHMMStateType
from src.antismash_parsing.genomic_context import ModuleGenomicContextFeature, ModuleGenomicContext
from math import log, e


def get_edge_weights(hmm,
                     hmm_scoring_cfg: HMMScoringConfig) -> Dict[Tuple[int, int], float]:  # hmm: DetailedHMM
    # for brevity
    ET = DetailedHMMEdgeType
    ST = DetailedHMMStateType
    MGF = ModuleGenomicContextFeature
    edge_weights_cfg = hmm_scoring_cfg.edge_weight_parameters

    # weights of the following edge types are determined based on the other edge types
    # these are the edge types corresponding to "natural" flow
    dependent_edge_types = {edge_type for edge_type in ET
                            if edge_type not in edge_weights_cfg}

    # q: assert that for each vertex among all transitions at most one has a dependent edge type
    # for debug
    assert all(
        len([edge for edge in hmm.transitions[u].values()
             if edge.edge_type in dependent_edge_types]) <= 1
        for u in range(len(hmm.states))
    )

    def get_edge_weight(_edge_info: DetailedHMMEdge):
        _edge_type = _edge_info.edge_type
        _edge_context = tuple(feature for feature in _edge_info.genomic_context
                             if feature in hmm_scoring_cfg.relevant_genomic_features[_edge_type])
        assert _edge_type in edge_weights_cfg and _edge_context in edge_weights_cfg[_edge_type], \
            f'Edge context {_edge_context} for {_edge_type} not found in edge weights config'
        return edge_weights_cfg[_edge_type][_edge_context]

    edge_weights = {}
    for u in range(len(hmm.states)):
        for v, edge_info in hmm.transitions[u].items():
            edge_type = edge_info.edge_type
            if edge_type in dependent_edge_types:
                pass
            elif edge_type in (ET.SKIP_MODULE_AT_START, ET.SKIP_MODULE_AT_END):
                    edge_weights[(u, v)] = 0.0  # penalties for skipping modules are added later
            else:
                edge_weights[(u, v)] = get_edge_weight(edge_info)

    # Accumulate skip-module penalties in the first edge of each skip-module path at start.
    for u in range(len(hmm.states)):
        if any(edge.edge_type == ET.START_MATCHING   # start of fragment skip at start
               for edge in hmm.transitions[u].values()):
            v = next((v for v in hmm.transitions[u]
                      if hmm.transitions[u][v].edge_type in (ET.START_SKIP_MODULES_AT_START,
                                                             ET.SKIP_MODULE_AT_START)),
                     None)
            if v is None:  # no more skip modules at start possible
                continue
            skip_path = [u, v]
            while len(hmm.transitions[v]) == 1:
                v = next(iter(hmm.transitions[v].keys()))
                skip_path.append(v)

            total_weight = sum(get_edge_weight(hmm.transitions[w1][w2])
                               for w1, w2 in pairwise(skip_path))
            # Move the total weight to the first edge in the path
            edge_weights[(skip_path[0], skip_path[1])] = total_weight

    # Accumulate skip-module penalties in the first edge of each skip-module path at end.
    for u in range(len(hmm.states)):
        v = next((v
                 for v, edge_info in hmm.transitions[u].items()
                 if edge_info.edge_type in (ET.END_MATCHING,
                                            ET.SKIP_MODULE_END_MATCHING)),
                 None)  # end of fragment skip at end):
        if v is None:
            continue

        w = next(w
                 for w, edge_info in hmm.transitions[v].items()
                 if hmm.states[w].state_type == ST.END_INSERTING_AT_END)

        skip_path = [w]
        while hmm.states[skip_path[-1]].state_type != ST.FINAL:
            assert len(hmm.transitions[skip_path[-1]]) == 1,\
                f'Expected only one transition from {hmm.states[skip_path[-1]].state_type}({skip_path[-1]})'
            skip_path.append(next(iter(hmm.transitions[skip_path[-1]].keys())))

        total_weight = sum(get_edge_weight(hmm.transitions[w1][w2])
                           for w1, w2 in pairwise(skip_path))
        # Move the total weight to the first edge in the path
        edge_weights[(u, v)] += total_weight


    # determine edge weights for the dependent edge types
    for u in range(len(hmm.states)):
        for v, edge_info in hmm.transitions[u].items():
            if edge_info.edge_type not in dependent_edge_types:
                continue

            sum_other_edges = sum(e ** edge_weights[(u, w)]
                                  for w in hmm.transitions[u]
                                  if w != v)
            try:
                edge_weights[(u, v)] = log(1 - sum_other_edges)
            except ValueError:
                raise ValueError(f'Invalid log probability for edge ({u}, {v}): '
                                 f'sum of other edges is {sum_other_edges}, '
                                 f'edge weights: {[edge_weights[(u, w)] for w in hmm.transitions[u] if w != v]}')

    # Assert that all edge weights are valid log probabilities
    eps = 1e-10  # small value to avoid numerical issues
    for u in range(len(hmm.states)):
        if len(hmm.transitions[u]) == 0:
            continue
        assert (all(edge_weights[(u, w)] <= 0.0 for w in hmm.transitions[u])
                and
                abs(1 - sum(e ** edge_weights[(u, w)] for w in hmm.transitions[u])) < eps), \
            f'Not all edge weights for vertex {u} are valid log probabilities: ' \
            f'{[edge_weights[(u, w)] for w in hmm.transitions[u]]}'

    return edge_weights
