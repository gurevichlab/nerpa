from typing import Dict, Tuple, TYPE_CHECKING, Union

from src.data_types import Prob

if TYPE_CHECKING:
    from src.matching.detailed_hmm import DetailedHMM

from src.matching.hmm_auxiliary_types import DetailedHMMEdgeType
from src.antismash_parsing.genomic_context import ModuleGenomicContextFeature
from src.matching.hmm_config import EdgeWeightsParams
from math import log


def get_edge_weights(hmm,  # type: DetailedHMM
                     probs: EdgeWeightsParams,
                     min_allowed_prob: float = 0) \
        -> Dict[Tuple[int, int], float]:
    # for brevity
    ET = DetailedHMMEdgeType
    MGF = ModuleGenomicContextFeature

    # weights of the following edge types are determined based on the other edge types
    # these are the edge types corresponding to "natural" flow
    dependent_edge_types = {
        ET.START_MATCHING,
        ET.MATCH,
        ET.NO_INSERTIONS,
        ET.END_INSERTING
    }

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
            edge_context = edge_info.genomic_context
            if edge_type in dependent_edge_types:
                continue

            edge_weights[(u, v)] = 0.1  # placeholder value

    '''
            def if_has_or_none(feature: MGF) -> Prob:
                return probs[edge_type][feature] if feature in edge_context \
                        else probs[edge_type][None]

            match edge_type:
                case (ET.START_INSERTING_AT_START | ET.INSERT_AT_START):
                    edge_weights[(u, v)] = if_has_or_none(MGF.PKS_UPSTREAM)
                case (ET.START_INSERTING | ET.INSERT):
                    pks_prob = probs[edge_type][MGF.PKS_DOWNSTREAM] \
                        if MGF.PKS_DOWNSTREAM in edge_context else 0
                    gene_end_prob = probs[edge_type][MGF.END_OF_GENE] \
                        if MGF.END_OF_GENE in edge_context else 0
                    default_prob = probs[edge_type][None]
                    edge_weights[(u, v)] = 1 - (1 - pks_prob) * (1 - gene_end_prob) * (1 - default_prob)
                case (ET.START_SKIP_MODULES_AT_START
                      | ET.SKIP_MODULE_AT_START
                      | ET.SKIP_MODULE
                      | ET.SKIP_MODULE_AT_END):
                    edge_weights[(u, v)] = if_has_or_none(MGF.ONLY_A_DOMAIN)

                case _:
                        edge_weights[(u, v)] = probs[edge_type][None]
    '''


    # determine edge weights for dependent edge types
    for u in range(len(hmm.states)):
        for v, edge_info in hmm.transitions[u].items():
            if edge_info.edge_type not in dependent_edge_types:
                continue

            sum_other_edges = sum(edge_weights[(u, w)]
                                  for w in hmm.transitions[u]
                                  if w != v)
            edge_weights[(u, v)] = 1 - sum_other_edges

    for u in range(len(hmm.states)):
        for v, edge_info in hmm.transitions[u].items():
            edge_weights[(u, v)] = log(edge_weights[(u, v)])
    return edge_weights
