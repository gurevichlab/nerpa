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
from src.training.training_types import MatchWithBGCNRP
from src.training.filter_edge_data import SingleFeatureContext, get_filtered_edge_data
from itertools import pairwise
from dataclasses import dataclass
from pathlib import Path
from collections import defaultdict
import yaml


@dataclass
class DataForTraining:
    edge_choices_cnts: Dict[DetailedHMMEdgeType, Dict[Tuple[SingleFeatureContext, ...], Tuple[int, int]]]
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
    #path = [state_idx for state_idx, _ in path_with_emissions]
    #detailed_hmm.draw(Path(f"{detailed_hmm.bgc_variant.genome_id}.png"),
    #                  highlight_path=path)
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
            start_state_idx = u - 1
            bgc_module = detailed_hmm.bgc_variant.modules[detailed_hmm.state_idx_to_module_idx[start_state_idx]]
            module_with_emission = (bgc_module, emission)
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

    edge_choices_cnts = get_filtered_edge_data(edge_choices)

    return DataForTraining(edge_choices_cnts=edge_choices_cnts,
                           match_emissions=match_emissions)
