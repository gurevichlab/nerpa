from typing import Dict, Tuple, TYPE_CHECKING, Union

from src.data_types import Prob, LogProb
from src.matching.hmm_config import HMMScoringConfig

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


    # determine edge weights for the dependent edge types
    for u in range(len(hmm.states)):
        for v, edge_info in hmm.transitions[u].items():
            if edge_info.edge_type not in dependent_edge_types:
                continue

            sum_other_edges = sum(e ** edge_weights[(u, w)]
                                  for w in hmm.transitions[u]
                                  if w != v)
            edge_weights[(u, v)] = log(1 - sum_other_edges)

    return edge_weights
