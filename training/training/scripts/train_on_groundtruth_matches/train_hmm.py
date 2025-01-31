from typing import (
    Dict,
    List,
    Tuple,
    Optional,
    Union
)
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
from src.monomer_names_helper import NRP_Monomer
from src.matching.matcher_viterbi_detailed_hmm import DetailedHMM
from src.write_results import write_yaml
from src.data_types import BGC_Module
from auxilary_types import MatchWithBGCNRP
from itertools import pairwise
from dataclasses import dataclass
from pathlib import Path
from collections import defaultdict
import yaml


@dataclass
class DataForTraining:
    edge_choices: List[Tuple[DetailedHMMEdgeType, GenomicContext, bool]]  # (edge_type, genomic_context, chosen)
    match_emissions:  List[Tuple[BGC_Module, NRP_Monomer]]  # (module, monomer)
    #insert_emissions: List[Tuple[BGC_Module, NRP_Monomer]]
    #insert_at_start_emissions: List[Tuple[BGC_Module, NRP_Monomer]]

@dataclass
class PathTurnInfo:
    chosen_edge_key: EdgeKey
    chosen_edge_info: Tuple[DetailedHMMEdgeType, GenomicContext]
    other_edges_info: List[Tuple[DetailedHMMEdgeType, GenomicContext]]
    bgc_module_with_emission: Optional[Tuple[BGC_Module, NRP_Monomer]] = None


def get_turns_info(detailed_hmm: DetailedHMM,
                   path_with_emissions: List[Tuple[int, Optional[NRP_Monomer]]]) \
        -> List[PathTurnInfo]:
    """
    At each turn of the path a choice between outcoming edges is made.
    This function returns the information about outcoming edges at each turn of the path,
    and indications whether the edge was chosen or not
    """
    num_insertions = 0
    turns_info = []
    for (u, emission), (v, _) in pairwise(path_with_emissions):
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

        if detailed_hmm.states[u].state_type == DetailedHMMStateType.MATCH:
            bgc_module = detailed_hmm.bgc_variant.modules[detailed_hmm.state_idx_to_module_idx[u]]
            emitted_monomer = path_with_emissions[u][1]
            module_with_emission = (bgc_module, emitted_monomer)
        else:
            module_with_emission = None

        turns_info.append(PathTurnInfo(chosen_edge_key,
                                       chosen_edge_info,
                                       other_edges_info,
                                       module_with_emission))

    return turns_info


def extract_data_for_training(matches_with_bgcs_nrps: List[MatchWithBGCNRP]) -> DataForTraining:
    turns_info = []
    for match, bgc_variant, nrp_variant in matches_with_bgcs_nrps:
        print(f"Processing match {match.nrp_variant_info.nrp_id}")
        detailed_hmm = DetailedHMM.from_bgc_variant(bgc_variant)
        # detailed_hmm.draw(Path(f"bgc.png"))
        for i, alignment in enumerate(match.alignments):
            #with open(f"alignment.txt", "w") as f:
            #    f.write(show_alignment(alignment))
            path_with_emissions = detailed_hmm.alignment_to_path_with_emisions(alignment)
            #path = [state_idx for state_idx, _ in path_with_emissions]
            #detailed_hmm.draw(Path(f"optimal_path.png"),
            #                  highlight_path=path)
            turns_info.extend(get_turns_info(detailed_hmm, path_with_emissions))

    seen_edge_keys = set()
    edge_choices = []
    match_emissions = []
    for turn_info in turns_info:
        if turn_info.chosen_edge_key in seen_edge_keys:
            continue

        edge_choices.append((turn_info.chosen_edge_info[0],
                             turn_info.chosen_edge_info[1],
                             True))
        if turn_info.bgc_module_with_emission is not None:
            match_emissions.append(turn_info.bgc_module_with_emission)

        for edge_info in turn_info.other_edges_info:
            edge_choices.append((edge_info[0], edge_info[1], False))

    return DataForTraining(edge_choices=edge_choices,
                           match_emissions=match_emissions)


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


EDGE_TYPE_DEPENDENCIES = {
    # Insertion edge dependencies
    DetailedHMMEdgeType.START_INSERTING_AT_START: {
        ModuleLocFeature.START_OF_BGC,
        ModuleLocFeature.START_OF_FRAGMENT,
        ModuleLocFeature.START_OF_GENE,
        ModuleLocFeature.PKS_UPSTREAM_PREV_GENE,
        ModuleLocFeature.PKS_UPSTREAM_SAME_GENE
    },
    DetailedHMMEdgeType.INSERT_AT_START: set(),
    DetailedHMMEdgeType.START_INSERTING: {
        ModuleLocFeature.END_OF_BGC,
        ModuleLocFeature.END_OF_FRAGMENT,
        ModuleLocFeature.END_OF_GENE,
        ModuleLocFeature.PKS_DOWNSTREAM_NEXT_GENE,
        ModuleLocFeature.PKS_DOWNSTREAM_SAME_GENE
    },
    DetailedHMMEdgeType.INSERT: set(),
    DetailedHMMEdgeType.END_INSERTING: set(),

    # Skip edge dependencies (empty for now)
    DetailedHMMEdgeType.START_SKIP_MODULES_AT_START: set(),
    DetailedHMMEdgeType.START_SKIP_GENES_AT_START: set(),
    DetailedHMMEdgeType.START_SKIP_FRAGMENTS_AT_START: set(),
    DetailedHMMEdgeType.SKIP_FRAGMENT_AT_START: set(),
    DetailedHMMEdgeType.SKIP_MODULE: set(),
    DetailedHMMEdgeType.SKIP_GENE: set(),
    DetailedHMMEdgeType.SKIP_FRAGMENT: set(),
    DetailedHMMEdgeType.START_SKIPPING_AT_END: set(),
    DetailedHMMEdgeType.SKIP_FRAGMENT_AT_END: set(),

    # Iteration edge dependencies (empty for now)
    DetailedHMMEdgeType.ITERATE_MODULE: set(),
    DetailedHMMEdgeType.ITERATE_GENE: set()
}

def get_filtered_edge_data(edge_choices: List[Tuple[DetailedHMMEdgeType, GenomicContext, bool]]) \
        -> Dict[DetailedHMMEdgeType, Dict[SingleFeatureContext, Tuple[int, int]]]:
    filtered_data = {}
    for edge_type, context, chosen in edge_choices:
        if edge_type not in EDGE_TYPE_DEPENDENCIES:
            continue
        if context is None:
            context = ()
        context = filter_context(context)
        context = tuple(feature for feature in context if feature in EDGE_TYPE_DEPENDENCIES[edge_type])
        if edge_type not in filtered_data:
            filtered_data[edge_type] = {}
        if context not in filtered_data[edge_type]:
            filtered_data[edge_type][context] = (0, 0)
        filtered_data[edge_type][context] = (filtered_data[edge_type][context][0] + 1 - int(chosen),
                                                filtered_data[edge_type][context][1] + int(chosen))

    write_yaml(filtered_data, Path('filtered_data.yaml'))
    return filtered_data


def estimate_parameters(matches_with_bgcs_nrps: List[MatchWithBGCNRP]) \
    -> Dict[DetailedHMMEdgeType, Dict[SingleFeatureContext, float]]:
    data = extract_data_for_training(matches_with_bgcs_nrps)
    edge_data = get_filtered_edge_data(data.edge_choices)
    edge_params = infer_edge_parameters(edge_data)
    emission_params = infer_emission_parameters(data.match_emissions)
    write_hmm_params(edge_params, emission_params)
    return edge_params, emission_params



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
