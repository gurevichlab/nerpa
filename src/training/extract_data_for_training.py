from typing import (
    Dict,
    List,
    NamedTuple,
    Set,
    Tuple,
    Optional,
    Union
)
from src.data_types import BGC_ID
from src.matching.hmm_auxiliary_types import (
    DetailedHMMStateType,
    DetailedHMMEdgeType,
    GenomicContext,
    StateIdx,
)
from src.matching.hmm_scoring_helper import HMMHelper
from src.monomer_names_helper import NRP_Monomer
from src.matching.detailed_hmm import DetailedHMM
from src.data_types import BGC_Module
from src.training.training_types import (
    MatchWithBGCNRP,
    DataForTraining,
    MatchEmissionInfo,
    ExtendedEdgeKey,
    EdgeInfo,
    PathTurnInfo, ChoicesCnts,
)
from src.training.filter_edge_data import get_filtered_edge_choices
from itertools import pairwise
from dataclasses import dataclass
from pathlib import Path
from collections import defaultdict
import yaml



def get_turns_info_for_match(detailed_hmm: DetailedHMM,
                             path_with_emissions: List[Tuple[StateIdx, Optional[NRP_Monomer]]]) \
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
        edge = detailed_hmm.transitions[u][v]
        if edge.edge_type == DetailedHMMEdgeType.INSERT:
            num_insertions += 1
        else:
            num_insertions = 0
        chosen_edge_key = ExtendedEdgeKey(edge.edge_key, num_insertions)

        chosen_edge_info = EdgeInfo(edge_type=edge.edge_type,
                                    genomic_context=edge.genomic_context)

        other_edges_info = []
        for w, edge_info in detailed_hmm.transitions[u].items():
            if w == v:
                continue
            other_edges_info.append(EdgeInfo(edge_type=edge_info.edge_type,
                                             genomic_context=edge_info.genomic_context))

        if detailed_hmm.states[u].state_type == DetailedHMMStateType.MATCH:
            try:
                next_module_start_state_idx = next(w
                                                   for w in detailed_hmm.transitions[u]
                                                   if detailed_hmm.states[w].state_type == DetailedHMMStateType.MODULE_START)
                next_module_idx = detailed_hmm.state_idx_to_module_idx[next_module_start_state_idx]
                module_idx = next_module_idx - 1
            except StopIteration:  # last module
                module_idx = len(detailed_hmm.bgc_variant.modules) - 1

            bgc_module = detailed_hmm.bgc_variant.modules[module_idx]
            module_with_emission = MatchEmissionInfo(bgc_module=bgc_module,
                                                     nrp_monomer=emission)
        else:
            module_with_emission = None

        turns_info.append(PathTurnInfo(chosen_edge_key,
                                       chosen_edge_info,
                                       other_edges_info,
                                       module_with_emission))

    return turns_info


def get_turns_info(matches_with_bgcs_nrps: List[MatchWithBGCNRP],
                   hmm_helper: HMMHelper) -> List[PathTurnInfo]:
    turns_info = []

    edge_keys_per_bgc: Dict[BGC_ID, Set[ExtendedEdgeKey]] = defaultdict(set)
    for match, bgc_variant, nrp_variant in matches_with_bgcs_nrps:
        print(f"Processing match {match.nrp_variant_id.nrp_id}")
        bgc_id = bgc_variant.bgc_variant_id.bgc_id
        detailed_hmm = DetailedHMM.from_bgc_variant(bgc_variant, hmm_helper)
        # detailed_hmm.draw(Path(f"bgc.png"))
        for i, alignment in enumerate(match.alignments):
            # with open(f"alignment.txt", "w") as f:
            #    f.write(show_alignment(alignment))
            path_with_emissions = detailed_hmm.alignment_to_path_with_emisions(alignment)
            # path = [state_idx for state_idx, _ in path_with_emissions]
            # detailed_hmm.draw(Path(f"optimal_path.png"),
            #                  highlight_path=path)
            new_turns = get_turns_info_for_match(detailed_hmm, path_with_emissions)
            for turn in new_turns:
                if turn.chosen_edge_key not in edge_keys_per_bgc[bgc_id]:
                    edge_keys_per_bgc[bgc_id].add(turn.chosen_edge_key)
                    turns_info.append(turn)

    return turns_info


def get_edge_choices(turns_info: List[PathTurnInfo]) -> List[Tuple[EdgeInfo, bool]]:
    edge_choices = []
    for turn_info in turns_info:
        edge_choices.append((turn_info.chosen_edge_info, True))
        edge_choices.extend((edge_info, False)
                            for edge_info in turn_info.other_edges_info)
    return edge_choices


def get_edge_choices_cnts(edge_choices: List[Tuple[EdgeInfo, bool]]) \
        -> Dict[DetailedHMMEdgeType, Dict[Optional[GenomicContext], ChoicesCnts]]:
    edge_choices_cnts = defaultdict(lambda: defaultdict(lambda: ChoicesCnts(0, 0)))
    for (edge_type, genomic_context), chosen in edge_choices:
        old_cnts = edge_choices_cnts[edge_type][genomic_context]
        edge_choices_cnts[edge_type][genomic_context] = ChoicesCnts(CHOSEN=old_cnts.CHOSEN + int(chosen),
                                                                    NOT_CHOSEN=old_cnts.NOT_CHOSEN + 1 - int(chosen))
    return edge_choices_cnts


def extract_data_for_training(matches_with_bgcs_nrps: List[MatchWithBGCNRP],
                              hmm_helper: HMMHelper) -> DataForTraining:
    turns_info = get_turns_info(matches_with_bgcs_nrps, hmm_helper)

    match_emissions = [turn_info.emission_info for turn_info in turns_info]
    edge_choices = get_edge_choices(turns_info)
    filtered_edge_choices = get_filtered_edge_choices(edge_choices)
    filtered_edge_choices_cnts = get_edge_choices_cnts(filtered_edge_choices)

    return DataForTraining(edge_choices_cnts=filtered_edge_choices_cnts,
                           match_emissions=match_emissions)
