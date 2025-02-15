from typing import Dict, Tuple, TYPE_CHECKING, Union

if TYPE_CHECKING:
    from src.matching.detailed_hmm import DetailedHMM

from src.matching.hmm_auxiliary_types import DetailedHMMEdgeType
from src.antismash_parsing.location_features import (
    ModuleLocFeature,
    GeneLocFeature,
    BGC_Fragment_Loc_Feature,
)
from src.matching.hmm_config import EdgeWeightsParams
from math import log


def get_edge_weights(hmm,  # type: DetailedHMM
                     probs: EdgeWeightsParams) \
        -> Dict[Tuple[int, int], float]:
    # weights of the following edge types are determined based on the other edge types
    # these are the edge types corresponding to "natural" flow
    dependent_edge_types = {
        DetailedHMMEdgeType.START_MATCHING,
        DetailedHMMEdgeType.MATCH,
        DetailedHMMEdgeType.NO_INSERTIONS,
        DetailedHMMEdgeType.END_INSERTING
    }

    # TODO: introduce other features, e.g. single module in gene
    genomic_location_features = {
        ModuleLocFeature.START_OF_BGC,
        ModuleLocFeature.END_OF_BGC,
        ModuleLocFeature.START_OF_FRAGMENT,
        ModuleLocFeature.END_OF_FRAGMENT,
        ModuleLocFeature.START_OF_GENE,
        ModuleLocFeature.END_OF_GENE,

        GeneLocFeature.START_OF_BGC,
        GeneLocFeature.END_OF_BGC,
        GeneLocFeature.START_OF_FRAGMENT,
        GeneLocFeature.END_OF_FRAGMENT,

        BGC_Fragment_Loc_Feature.START_OF_BGC,
        BGC_Fragment_Loc_Feature.END_OF_BGC
    }

    pks_context_features = {
        ModuleLocFeature.PKS_UPSTREAM_PREV_GENE,
        ModuleLocFeature.PKS_DOWNSTREAM_NEXT_GENE,
        ModuleLocFeature.PKS_UPSTREAM_SAME_GENE,
        ModuleLocFeature.PKS_DOWNSTREAM_SAME_GENE,

        GeneLocFeature.PKS_UPSTREAM,
        GeneLocFeature.PKS_DOWNSTREAM,

        BGC_Fragment_Loc_Feature.PKS_UPSTREAM,
        BGC_Fragment_Loc_Feature.PKS_DOWNSTREAM
    }

    # determine edge weights for non-dependent edge types
    edge_weights = {}
    for u in range(len(hmm.states)):
        for v, edge_info in hmm.adj_list[u].items():
            if edge_info.edge_type in dependent_edge_types:
                continue

            if edge_info.genomic_context is None:
                genomic_loc_prob = 0
            else:
                genomic_loc_prob = max((probs[edge_info.edge_type][loc_feature]
                                    for loc_feature in genomic_location_features
                                    if loc_feature in edge_info.genomic_context
                                    and loc_feature in probs[edge_info.edge_type]),
                                   default=0)

            if edge_info.genomic_context is None:
                pks_context_prob = 0
            else:
                pks_context_prob = max((probs[edge_info.edge_type][loc_feature]
                                        for loc_feature in pks_context_features
                                        if loc_feature in edge_info.genomic_context
                                        and loc_feature in probs[edge_info.edge_type]),
                                       default=0)
            base_prob = probs[edge_info.edge_type].get(None, 1e-3)
            edge_weights[(u, v)] = 1 - (1 - genomic_loc_prob) * (1 - pks_context_prob) * (1 - base_prob)

    # determine edge weights for dependent edge types
    for u in range(len(hmm.states)):
        for v, edge_info in hmm.adj_list[u].items():
            if edge_info.edge_type not in dependent_edge_types:
                continue

            sum_other_edges = sum(edge_weights[(u, w)]
                                  for w in hmm.adj_list[u]
                                  if w != v)
            edge_weights[(u, v)] = 1 - sum_other_edges

    for u in range(len(hmm.states)):
        for v, edge_info in hmm.adj_list[u].items():
            edge_weights[(u, v)] = log(edge_weights[(u, v)])
    return edge_weights
