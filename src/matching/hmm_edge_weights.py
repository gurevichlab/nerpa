from itertools import pairwise
from typing import Dict, Tuple, TYPE_CHECKING, Union

from src.data_types import Prob, LogProb
from src.matching.hmm_scoring_config import HMMScoringConfig

if TYPE_CHECKING:
    from src.matching.detailed_hmm import DetailedHMM

from src.matching.hmm_auxiliary_types import DetailedHMMEdgeType
from src.antismash_parsing.genomic_context import ModuleGenomicContextFeature, ModuleGenomicContext
from math import log, e


def get_edge_weights(hmm,
                     hmm_scoring_cfg: HMMScoringConfig) -> Dict[Tuple[int, int], float]:  # hmm: DetailedHMM
    # for brevity
    ET = DetailedHMMEdgeType
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

    edge_weights = {}
    for u in range(len(hmm.states)):
        for v, edge_info in hmm.transitions[u].items():
            edge_type = edge_info.edge_type
            if edge_type in dependent_edge_types:
                continue

            edge_context = tuple(feature for feature in edge_info.genomic_context
                                 if feature in hmm_scoring_cfg.relevant_genomic_features[edge_type])
            assert edge_type in edge_weights_cfg and edge_context in edge_weights_cfg[edge_type], \
                f'Edge context {edge_context} for {edge_type} not found in edge weights config'
            edge_weights[(u, v)] = edge_weights_cfg[edge_type][edge_context]

    # Move total weights of non-branching paths to the first edge in the path
    # this is done to ensure that weights of transitions from each state correspond to a probability distribution
    for u in range(len(hmm.states)):
        if len(hmm.transitions[u]) > 1:
            for v in hmm.transitions[u]:
                if len(hmm.transitions[v]) != 1:
                    continue
                non_branching_path = [u, v]
                w = v
                while len(hmm.transitions[w]) == 1:
                    edge_end = next(iter(hmm.transitions[w].keys()))
                    non_branching_path.append(edge_end)
                    w = edge_end

                total_weight = sum(edge_weights[edge]
                                   for edge in pairwise(non_branching_path))
                # Move the total weight to the first edge in the path
                edge_weights[(u, v)] = total_weight
                for not_fst_edge in pairwise(non_branching_path)[1:]:
                    edge_weights[not_fst_edge] = 0.0



    # determine edge weights for the dependent edge types
    for u in range(len(hmm.states)):
        for v, edge_info in hmm.transitions[u].items():
            if edge_info.edge_type not in dependent_edge_types:
                continue

            sum_other_edges = sum(e ** edge_weights[(u, w)]
                                  for w in hmm.transitions[u]
                                  if w != v)
            edge_weights[(u, v)] = log(1 - sum_other_edges)

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
