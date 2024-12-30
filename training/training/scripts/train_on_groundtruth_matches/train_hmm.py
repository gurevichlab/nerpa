from typing import Dict, List, Tuple, Union
from src.antismash_parsing.location_features import ModuleLocFeature, GeneLocFeature, BGC_Fragment_Loc_Feature
from src.matching.matcher_viterbi_types import (
    HMM,
    DetailedHMMStateType,
    DetailedHMMState,
    DetailedHMMEdgeType,
    DetailedHMMEdge,
    GenomicContext,
    EdgeKey
)
from src.matching.matching_types_alignment import show_alignment
from src.matching.matcher_viterbi_detailed_hmm import DetailedHMM
from auxilary_types import MatchWithBGCNRP
from itertools import pairwise
from dataclasses import dataclass
from pathlib import Path
from collections import defaultdict


@dataclass
class DataForTraining:
    edge_choices: List[Tuple[DetailedHMMEdgeType, GenomicContext, bool]]  # (edge_type, genomic_context, chosen)
    #emissions  # TBWR

@dataclass
class PathTurnInfo:
    chosen_edge_key: EdgeKey
    chosen_edge_info: Tuple[DetailedHMMEdgeType, GenomicContext]
    other_edges_info: List[Tuple[DetailedHMMEdgeType, GenomicContext]]


def get_turns_info(detailed_hmm: DetailedHMM, path: List[int]) \
        -> List[PathTurnInfo]:
    """
    At each turn of the path a choice between outcoming edges is made.
    This function returns the information about outcoming edges at each turn of the path,
    and indications whether the edge was chosen or not
    """
    num_insertions = 0
    turns_info = []
    for u, v in pairwise(path):
        edge = detailed_hmm.adj_list[u][v]
        if edge.edge_type == DetailedHMMEdgeType.INSERT:
            chosen_edge_key = edge.edge_key + (num_insertions,)
            num_insertions += 1
        else:
            chosen_edge_key = edge.edge_key
            num_insertions = 0

        chosen_edge_info = (edge.edge_type, edge.genomic_context)

        other_edges_info = []
        for w, edge_info in detailed_hmm.adj_list[u].items():
            if w == v:
                continue
            other_edges_info.append((edge_info.edge_type, edge_info.genomic_context))

        turns_info.append(PathTurnInfo(chosen_edge_key, chosen_edge_info, other_edges_info))

    return turns_info


def extract_data_for_training(matches_with_bgcs_nrps: List[MatchWithBGCNRP]) -> DataForTraining:
    turns_info = []
    for match, bgc_variant, nrp_variant in matches_with_bgcs_nrps:
        detailed_hmm = DetailedHMM.from_bgc_variant(bgc_variant)
        detailed_hmm.draw(Path(f"bgc_{bgc_variant.genome_id}.png"))
        for i, alignment in enumerate(match.alignments):
            with open(f"bgc_{bgc_variant.genome_id}_alignment_{i}.txt", "w") as f:
                f.write(show_alignment(alignment))
            path_with_emissions = detailed_hmm.alignment_to_path_with_emisions(alignment)
            path = [state_idx for state_idx, _ in path_with_emissions]
            detailed_hmm.draw(Path(f"bgc_{bgc_variant.genome_id}_alignment_{i}.png"),
                              highlight_path=path)
            turns_info.extend(get_turns_info(detailed_hmm, path))

    seen_edge_keys = set()
    edge_choices = []
    for turn_info in turns_info:
        if turn_info.chosen_edge_key in seen_edge_keys:
            continue

        edge_choices.append((turn_info.chosen_edge_info[0],
                             turn_info.chosen_edge_info[1],
                             True))
        for edge_info in turn_info.other_edges_info:
            edge_choices.append((edge_info[0], edge_info[1], False))

    return DataForTraining(edge_choices=edge_choices)


def filter_context(genomic_context: GenomicContext) -> GenomicContext:
    # if ModuleLocFeature.START_OF_BGC is present, then MODULE_LOC_FEATURE.START_OF_FRAGMENT is automatically present
    # I want to keep only ModuleLocFeature.START_OF_BGC.
    # The same for other pairs of features
    if ModuleLocFeature.START_OF_FRAGMENT in genomic_context:
        genomic_context = tuple(loc for loc in genomic_context if loc != ModuleLocFeature.START_OF_GENE)
    if ModuleLocFeature.END_OF_FRAGMENT in genomic_context:
        genomic_context = tuple(loc for loc in genomic_context if loc != ModuleLocFeature.END_OF_GENE)
    if ModuleLocFeature.START_OF_BGC in genomic_context:
        genomic_context = tuple(loc for loc in genomic_context if loc != ModuleLocFeature.START_OF_FRAGMENT)
    if ModuleLocFeature.END_OF_BGC in genomic_context:
        genomic_context = tuple(loc for loc in genomic_context if loc != ModuleLocFeature.END_OF_FRAGMENT)
    return genomic_context


# parameters are estimated for single features only
# for multiple features, parameters are computed based on the single features
SingleFeatureContext = Union[ModuleLocFeature, GeneLocFeature, BGC_Fragment_Loc_Feature, None]

def estimate_parameters(matches_with_bgcs_nrps: List[MatchWithBGCNRP]) \
    -> Dict[DetailedHMMEdgeType, Dict[SingleFeatureContext, float]]:
    data = extract_data_for_training(matches_with_bgcs_nrps)
    cnts = defaultdict(lambda: defaultdict(lambda: (0,0)))
    for edge_type, genomic_context, chosen in data.edge_choices:
        if not genomic_context:
            cnts[edge_type][None] = (cnts[edge_type][genomic_context][0] + 1 - int(chosen),
                                     cnts[edge_type][genomic_context][1] + int(chosen))
            continue
        genomic_context = filter_context(genomic_context)
        for loc_feature in genomic_context:
            cnts[edge_type][loc_feature] = (cnts[edge_type][loc_feature][0] + 1 - int(chosen),
                                            cnts[edge_type][loc_feature][1] + int(chosen))

    default_prob = {}
    for edge_type in cnts:
        total_not_chosen = sum(cnts[edge_type][loc_feature][0] for loc_feature in cnts[edge_type])
        total_chosen = sum(cnts[edge_type][loc_feature][1] for loc_feature in cnts[edge_type])

        default_prob[edge_type] = total_chosen / (total_chosen + total_not_chosen)

    probs = defaultdict(lambda: defaultdict(float))
    for edge_type in cnts:
        for loc_feature in cnts[edge_type]:
            not_chosen, chosen = cnts[edge_type][loc_feature]
            if not_chosen + chosen < 10:
                print(f"Warning: not enough data for edge type {edge_type} and loc feature {loc_feature}")

            probs[edge_type][loc_feature] = chosen / (chosen + not_chosen)
    return probs


def get_edge_weights(hmm: DetailedHMM,
                     probs: Dict[DetailedHMMEdgeType, Dict[SingleFeatureContext, float]]) \
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

            genomic_loc_prob = max((probs[edge_info.edge_type][loc_feature]
                                    for loc_feature in genomic_location_features
                                    if loc_feature in edge_info.genomic_context),
                                   default=0)
            pks_context_prob = max((probs[edge_info.edge_type][loc_feature]
                                    for loc_feature in pks_context_features
                                    if loc_feature in edge_info.genomic_context),
                                   default=0)
            base_prob = probs[edge_info.edge_type][None]
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

    return edge_weights


'''
mu_p = 0.1  # prior mean TODO: base it on default_prob[edge_type] or smth instead of hardcoding
k = 1  # concentration parameter (k = alpha + beta for Beta distribution)
prior_params = (mu_p, 1 - mu_p)  # Beta distribution parameters
posterior_params = (prior_params[0] + chosen, prior_params[1] + not_chosen)
posterior_mean = posterior_params[0] / sum(posterior_params)

probs[edge_type][loc_feature] = posterior_mean
'''
