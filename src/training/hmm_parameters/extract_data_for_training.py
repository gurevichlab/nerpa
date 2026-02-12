from pathlib import Path
from typing import (
    Dict,
    List,
    Tuple,
    Optional
)

import pandas as pd

from src.antismash_parsing.antismash_parser_types import BGC_ID
from src.matching.alignment_type import light_alignment_to_str
from src.hmm.hmm_auxiliary_types import (
    DetailedHMMStateType,
    DetailedHMMEdgeType,
    GenomicContext,
    StateIdx,
)
from src.hmm.hmm_constructor.hmm_constructor_state_edge_context_relations import MATCHING_STATE_TYPES, RELEVANT_GENOMIC_CONTEXT
from src.hmm.hmm_scoring_helper import HMMHelper
from src.hmm.detailed_hmm import DetailedHMM
from src.rban_parsing.nrp_variant_types import NRP_Variant
from src.rban_parsing.rban_monomer import rBAN_Monomer
from src.testing.simplified_alignment import simplified_alignment_to_light_alignments, simplified_alignment_to_str
from src.training.hmm_parameters.training_types import (
    MatchWithBGCNRP,
    DataForTraining,
    EmissionInfo,
    EdgeInfo,
    PathTurnInfo,
    ChoicesCnts,
    EmissionKey,
    EdgeKey,
    EdgeChoices_df,
    EdgeChoicesSchema,
    make_edge_choices_df,
)
from itertools import pairwise
from collections import defaultdict

PathWithEmissions = List[Tuple[StateIdx, Optional[rBAN_Monomer]]]

def get_emissions_for_match(detailed_hmm: DetailedHMM,
                            path_with_emissions: PathWithEmissions) \
        -> List[EmissionInfo]:
    """
    Get the emissions in the path
    """
    emissions = []
    match_state_idx_to_module_idx = {state_idx: module_idx
                                     for module_idx, state_idx in enumerate(detailed_hmm._module_idx_to_match_state_idx)}

    for state_idx, emission in path_with_emissions:
        state_type = detailed_hmm.states[state_idx].state_type
        if emission is not None:
            if state_type in MATCHING_STATE_TYPES:
                module_idx = match_state_idx_to_module_idx[state_idx]
                bgc_module = detailed_hmm.bgc_variant.modules[module_idx]
            else:
                bgc_module = None
            emissions.append(EmissionInfo(bgc_id=detailed_hmm.bgc_variant.bgc_variant_id.bgc_id,
                                          bgc_module=bgc_module,
                                          nrp_monomer=emission,
                                          state_type=state_type))
    return emissions


def get_emissions(hmms_with_paths_with_emissions: List[Tuple[DetailedHMM, PathWithEmissions]]) \
        -> List[EmissionInfo]:
    #print('get_emissions: ', [path for hmm, path in hmms_with_paths_with_emissions])
    emissions = []
    seen_emissions = set()
    for detailed_hmm, path_with_emissions in hmms_with_paths_with_emissions:
        new_emissions = get_emissions_for_match(detailed_hmm, path_with_emissions)
        for emission in new_emissions:
            emission_key = EmissionKey.from_emission_info(emission)
            if emission_key not in seen_emissions:
                seen_emissions.add(emission_key)
                emissions.append(emission)

    #print('get_emissions result: ', emissions)
    return emissions


def get_turns_info_for_match(detailed_hmm: DetailedHMM,
                             path_with_emissions: PathWithEmissions) \
        -> List[PathTurnInfo]:
    """
    At each turn of the path a choice between outcoming edges is made.
    This function returns the information about outcoming edges at each turn of the path,
    and indications whether the edge was chosen or not
    """
    #path = [state_idx for state_idx, _ in path_with_emissions]
    #detailed_hmm.draw(Path(f"{detailed_hmm.bgc_variant.genome_id}.png"),
    #                  highlight_path=path)
    ET = DetailedHMMEdgeType
    bgc_id = detailed_hmm.bgc_variant.bgc_variant_id.bgc_id

    num_insertions = 0
    turns_info = []
    for (u, emission), (v, _) in pairwise(path_with_emissions):
        state = detailed_hmm.states[u]
        module = detailed_hmm.bgc_variant.modules[state.related_module_idx] \
            if state.related_module_idx is not None else None
        edge = detailed_hmm.transitions[u][v]
        num_insertions = num_insertions + 1 \
             if edge.edge_type in (ET.INSERT, ET.INSERT_AT_START) else 0

        chosen_edge_key = EdgeKey(u, v, num_insertions)
        chosen_edge_info = EdgeInfo(edge_type=edge.edge_type,
                                    edge_start_type=state.state_type,
                                    genomic_context=module.genomic_context)

        other_edges_info = []
        for w, edge_info in detailed_hmm.transitions[u].items():
            if w == v:  # skip the chosen edge
                continue
            other_edges_info.append(EdgeInfo(edge_type=edge_info.edge_type,
                                             edge_start_type=detailed_hmm.states[u].state_type,
                                             genomic_context=module.genomic_context))

        turns_info.append(PathTurnInfo(bgc_id,
                                       chosen_edge_key,
                                       chosen_edge_info,
                                       other_edges_info))

    return turns_info


def get_edge_choices_for_match(detailed_hmm: DetailedHMM,
                               nrp_variant: NRP_Variant,
                               path_with_emissions: PathWithEmissions) \
        -> EdgeChoices_df:
    """
    At each turn of the path a choice between outcoming edges is made.
    This function returns the information about outcoming edges at each turn of the path,
    and indications whether the edge was chosen or not
    """
    #path = [state_idx for state_idx, _ in path_with_emissions]
    #detailed_hmm.draw(Path(f"{detailed_hmm.bgc_variant.genome_id}.png"),
    #                  highlight_path=path)
    ET = DetailedHMMEdgeType
    ECS = EdgeChoicesSchema

    num_insertions = 0
    edge_choices = []
    for (u, emission), (v, _) in pairwise(path_with_emissions):
        state = detailed_hmm.states[u]
        module = detailed_hmm.bgc_variant.modules[state.related_module_idx] \
            if state.related_module_idx is not None else None
        edge = detailed_hmm.transitions[u][v]
        num_insertions = num_insertions + 1 \
            if edge.edge_type in (ET.INSERT, ET.INSERT_AT_START) else 0

        edge_choices.append({
            ECS.BGC_ID: detailed_hmm.bgc_variant.bgc_variant_id,
            ECS.NRP_ID: nrp_variant.nrp_variant_id,
            ECS.STATE_TYPE: state.state_type,
            ECS.EDGE_TYPE: edge.edge_type,
            ECS.NUM_INSERTIONS: num_insertions,
            ECS.EDGE_FROM: u,
            ECS.EDGE_TO: v,
            ECS.MODULE: module,
            ECS.EMISSION: emission,
        })



    return make_edge_choices_df(pd.DataFrame(edge_choices))


def get_edge_choices(hmms_with_paths_with_emissions: List[Tuple[DetailedHMM, NRP_Variant, PathWithEmissions]]) -> EdgeChoices_df:
    edge_choices = pd.DataFrame()
    # I used to filter duplicate edges per BGC, but now I want to keep all occurrences
    # edge_keys_per_bgc: Dict[BGC_ID, Set[ExtendedEdgeKey]] = defaultdict(set)

    for nrp_variant, detailed_hmm, path_with_emissions in hmms_with_paths_with_emissions:
        # bgc_id = detailed_hmm.bgc_variant.bgc_variant_id.bgc_id
        new_edge_choices = get_edge_choices_for_match(nrp_variant, detailed_hmm, path_with_emissions)
        edge_choices = pd.concat([edge_choices, new_edge_choices], ignore_index=True)
        # for turn in new_turns:
        #     if turn.chosen_edge_key not in edge_keys_per_bgc[bgc_id]:
        #         edge_keys_per_bgc[bgc_id].add(turn.chosen_edge_key)
        #         turns_info.append(turn)

    return make_edge_choices_df(edge_choices)


def get_hmms_with_paths_with_emissions(matches_with_bgcs_nrps: List[MatchWithBGCNRP],
                                       hmm_helper: HMMHelper,
                                       dir_for_hmm_figures: Optional[Path] = None) -> List[Tuple[DetailedHMM, NRP_Variant, PathWithEmissions]]:
    hmms_with_paths_with_emissions = []
    for match, bgc_variant, nrp_variant in matches_with_bgcs_nrps:
        try:
            detailed_hmm = DetailedHMM.from_bgc_variant(bgc_variant, hmm_helper)
        except Exception as e:
            print(f'WARNING: cannot construct DetailedHMM for match {match.nrp_id} ')
            raise
        if match.nrp_id in ('BGC0000395.0',   # Abu is trimmed as an isolated unknown monomer; however, it should be aligned
                            'BGC0000416.0', 'BGC0000416.1',  # gene XNC1_2229 is iterated although there're no signes of it
                            ):
            continue
        alignments_iter = iter(simplified_alignment_to_light_alignments(match.true_alignment,
                                                                        detailed_hmm.bgc_variant,
                                                                        nrp_variant))
        all_alignments = list(alignments_iter)
        if not all_alignments:
            print(f'WARNING: cannot reconstruct alignment for match {match.nrp_id} ')
            print('Original alignment:')
            print(simplified_alignment_to_str(match.true_alignment))
            print('NRP variant:')
            print(nrp_variant.to_str_compact())
            raise ValueError(f'Cannot reconstruct alignment for match {match.nrp_id} ')

        '''
        if (other_alignments := next(alignments_iter, None)) is not None:
            print(f'WARNING: multiple alignments for match {match.nrp_id} ')
            print('First alignment:')
            for al in alignments:
                print(light_alignment_to_str(al))
            print('Other alignment:')
            for al in other_alignments:
                print(light_alignment_to_str(al))
            raise ValueError(f'Ambiguous alignment reconstruction for match {match.nrp_id} '
                             f'after simplification to light alignments')
        '''
        if len(all_alignments) > 1:
            print(f'WARNING: multiple alignments for match {match.nrp_id} ')
            print('True alignment:')
            print(simplified_alignment_to_str(match.true_alignment))
            print('First reconstructed alignment:')
            for al in all_alignments[0]:
                print(light_alignment_to_str(al))
            print('Other reconstructed alignment:')
            for al in all_alignments[1]:
                print(light_alignment_to_str(al))
            print(f'WARNING: ambiguous alignment reconstruction for match {match.nrp_id} '
                  'picking the first one')
            # raise ValueError(f'Ambiguous alignment reconstruction for match {match.nrp_id} '
            #                  f'after simplification to light alignments')

        alignments = all_alignments[0]
        for i, alignment in enumerate(alignments):
            try:
                path_with_emissions = detailed_hmm.alignment_to_path_with_emisions(alignment,
                                                                                   alignment_type='light')
                #print('extract_data: ', path_with_emissions)
                if dir_for_hmm_figures is not None:
                    detailed_hmm.draw(dir_for_hmm_figures / Path(f'{detailed_hmm.bgc_variant.bgc_variant_id.bgc_id.antiSMASH_file}_path_{i}.png'),
                                      highlight_path=[state for state, _ in path_with_emissions])  # for debugging
            except Exception as e:
                print(f'WARNING: alignment {i} for match {match.nrp_id} '
                      f'could not be converted to path with emissions')
                raise e
            hmms_with_paths_with_emissions.append((detailed_hmm, nrp_variant, path_with_emissions))

    return hmms_with_paths_with_emissions


def turns_info_to_edge_choices(turns_info: List[PathTurnInfo]) -> List[Tuple[BGC_ID, EdgeInfo, bool]]:
    edge_choices = []
    for turn_info in turns_info:
        edge_choices.append((turn_info.bgc_id, turn_info.chosen_edge_info, True))
        edge_choices.extend((turn_info.bgc_id, edge_info, False)
                            for edge_info in turn_info.other_edges_info)
    return edge_choices


def get_edge_choices_cnts(edge_choices: List[Tuple[BGC_ID, EdgeInfo, bool]]) \
        -> Dict[DetailedHMMEdgeType, Dict[Optional[GenomicContext], ChoicesCnts]]:
    edge_choices_cnts = defaultdict(lambda: defaultdict(lambda: ChoicesCnts(0, 0)))
    for bgc_id, (edge_type, genomic_context), chosen in edge_choices:
        old_cnts = edge_choices_cnts[edge_type][genomic_context]
        edge_choices_cnts[edge_type][genomic_context] = ChoicesCnts(CHOSEN=old_cnts.CHOSEN + int(chosen),
                                                                    NOT_CHOSEN=old_cnts.NOT_CHOSEN + 1 - int(chosen))
    return edge_choices_cnts


def get_turns_info(hmms_with_paths_with_emissions: List[Tuple[DetailedHMM, PathWithEmissions]]) -> List[PathTurnInfo]:
    turns_info = []
    # I used to filter duplicate edges per BGC, but now I want to keep all occurrences
    # edge_keys_per_bgc: Dict[BGC_ID, Set[ExtendedEdgeKey]] = defaultdict(set)

    for detailed_hmm, path_with_emissions in hmms_with_paths_with_emissions:
        # bgc_id = detailed_hmm.bgc_variant.bgc_variant_id.bgc_id
        new_turns = get_turns_info_for_match(detailed_hmm, path_with_emissions)
        turns_info.extend(new_turns)
        # for turn in new_turns:
        #     if turn.chosen_edge_key not in edge_keys_per_bgc[bgc_id]:
        #         edge_keys_per_bgc[bgc_id].add(turn.chosen_edge_key)
        #         turns_info.append(turn)


    return turns_info


def get_turn_cnts_per_state_cntxt(turns_info: List[PathTurnInfo],
                                  filter_irrelevant_context: bool = True) \
        -> Dict[DetailedHMMStateType, Dict[GenomicContext, Dict[DetailedHMMEdgeType, int]]]:
    """
    Counts how many turns (choices between edges) are made for each edge type
     from each state type and genomic context.
    """
    # 1. group turns by state type and genomic context
    turns_per_state_cntxt: Dict[Tuple[DetailedHMMStateType, GenomicContext], List[DetailedHMMEdgeType]] = \
        defaultdict(list)
    for turn_info in turns_info:
        state_type = turn_info.chosen_edge_info.edge_start_type
        genomic_context = turn_info.chosen_edge_info.genomic_context

        turns_per_state_cntxt[(state_type, genomic_context)].append(turn_info.chosen_edge_info.edge_type)

    # 2. remove irrelevant context
    if filter_irrelevant_context:
        turns_per_state_cntxt_filtered: Dict[Tuple[DetailedHMMStateType, GenomicContext], List[DetailedHMMEdgeType]] = \
            defaultdict(list)
        for (state_type, genomic_context), edge_types in turns_per_state_cntxt.items():
            # for each edge type, keep only relevant genomic context features
            relevant_context = RELEVANT_GENOMIC_CONTEXT[state_type]
            filtered_context = tuple(feature
                                     for feature in genomic_context
                                     if feature in relevant_context)
            turns_per_state_cntxt_filtered[(state_type, filtered_context)].extend(edge_types)
    else:
        turns_per_state_cntxt_filtered = turns_per_state_cntxt

    # 3. count edge types for each (state type, genomic context)
    turns_per_state_cntxt_cnts: Dict[DetailedHMMStateType, Dict[GenomicContext, Dict[DetailedHMMEdgeType, int]]] = \
        defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for (state_type, genomic_context), edge_types in turns_per_state_cntxt_filtered.items():
        for edge_type in edge_types:
            turns_per_state_cntxt_cnts[state_type][genomic_context][edge_type] += 1

    return turns_per_state_cntxt_cnts


def extract_data_for_training(matches_with_bgcs_nrps: List[MatchWithBGCNRP],
                              hmm_helper: HMMHelper,
                              dir_for_hmm_figures: Optional[Path] = None) -> DataForTraining:
    hmms_nrps_paths = get_hmms_with_paths_with_emissions(matches_with_bgcs_nrps,
                                                                        hmm_helper,
                                                                        dir_for_hmm_figures)
    hmms_paths = [(hmm, path) for hmm, _, path in hmms_nrps_paths]
    edge_choices_df = get_edge_choices(hmms_nrps_paths)

    turns_info = get_turns_info(hmms_paths)
    turns_per_state_cntxt_cnts = get_turn_cnts_per_state_cntxt(turns_info)
    turns_per_state_cntxt_cnts_wo_filt = get_turn_cnts_per_state_cntxt(turns_info,
                                                                       filter_irrelevant_context=False)
    emissions = get_emissions(hmms_paths)

    '''
    edge_choices = turns_info_to_edge_choices(turns_info)
    filtered_edge_choices = get_filtered_edge_choices(edge_choices)
    filtered_edge_choices_cnts = get_edge_choices_cnts(filtered_edge_choices)

    chosen_edges_occurrences = defaultdict(lambda: defaultdict(list))
    for bgc_id, edge_info, chosen in filtered_edge_choices:
        if chosen:
            chosen_edges_occurrences[edge_info.edge_type][edge_info.genomic_context].append(bgc_id)
    '''

    return DataForTraining(edge_choices_per_state_cntxt=turns_per_state_cntxt_cnts,
                           emissions=emissions,
                           edge_choices_wo_filtering=turns_per_state_cntxt_cnts_wo_filt,
                           edge_choices_df=edge_choices_df)


