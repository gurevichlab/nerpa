from typing import (
    Dict,
    Iterable,
    List,
    NamedTuple
)
from auxilary_types import MatchWithBGCNRP
from dataclasses import dataclass
from src.antismash_parsing.location_features import (
    ModuleLocFeatures,
    ModuleLocFeature,
    GeneLocFeatures,
    GeneLocFeature,
    BGC_Fragment_Loc_Features,
    BGC_Fragment_Loc_Feature
)
from src.matching.matching_types_alignment_step import (
    AlignmentStep,
    AlignmentStep_BGC_Module_Info
)
from src.matching.matching_types_match import Match
from src.monomer_names_helper import MonomerResidue
from src.data_types import (
    BGC_Variant,
    GeneId,
    LogProb,
    NRP_Variant
)
from extract_insert_steps import InsertRunsInfo, extract_inserts_info
from extract_matches_skips_steps import MatchesSkipsSteps, extract_matches_skips_info


@dataclass
class StepsInfo:
    matches_and_skips: MatchesSkipsSteps
    insert_runs: InsertRunsInfo


def get_steps_info(matches_with_bgcs_nrps: List[MatchWithBGCNRP]) -> StepsInfo:
    matches_and_skips = MatchesSkipsSteps(module_matches=[], module_skips=[],
                                          gene_matches=[], gene_skips=[],
                                          bgc_fragment_matches=[], bgc_fragment_skips=[])
    insert_runs = InsertRunsInfo(inserts=[], inserts_at_start=[])

    # q: group matches by bgc_id
    bgc_id_to_bgc_variant = {(match_with_bgcs_nrps.bgc_variant.genome_id, match_with_bgcs_nrps.bgc_variant.bgc_idx): match_with_bgcs_nrps.bgc_variant
                              for match_with_bgcs_nrps in matches_with_bgcs_nrps}
    bgc_id_to_matches = {(match_with_bgcs_nrps.bgc_variant.genome_id, match_with_bgcs_nrps.bgc_variant.bgc_idx): []
                         for match_with_bgcs_nrps in matches_with_bgcs_nrps}
    for match_with_bgcs_nrps in matches_with_bgcs_nrps:
        bgc_id_to_matches[(match_with_bgcs_nrps.bgc_variant.genome_id, match_with_bgcs_nrps.bgc_variant.bgc_idx)].append(match_with_bgcs_nrps.match)

    for bgc_id, bgc_variant in bgc_id_to_bgc_variant.items():
        matches = bgc_id_to_matches[bgc_id]
        new_steps_info = get_matches_stats_for_bgc(matches, bgc_variant)
        matches_and_skips.join(new_steps_info.matches_and_skips, unique=False)
        insert_runs.join(new_steps_info.insert_runs, unique=False)
    return StepsInfo(matches_and_skips, insert_runs)


# TODO: in the future normalize counts for each NRP family to avoid overfitting
def get_matches_stats_for_bgc(matches: List[Match],
                              bgc_variant: BGC_Variant) -> StepsInfo:
    matches_and_skips = MatchesSkipsSteps(module_matches=[], module_skips=[],
                                          gene_matches=[], gene_skips=[],
                                          bgc_fragment_matches=[], bgc_fragment_skips=[])
    insert_runs = InsertRunsInfo(inserts=[], inserts_at_start=[])
    for match in matches:
        for alignment in match.alignments:
            new_matches_and_skips = extract_matches_skips_info(alignment, bgc_variant, match.nrp_variant_info.nrp_id)
            new_insert_runs = extract_inserts_info(alignment, bgc_variant, match.nrp_variant_info.nrp_id)
            matches_and_skips.join(new_matches_and_skips, unique=False)
            insert_runs.join(new_insert_runs, unique=False)
    return StepsInfo(matches_and_skips, insert_runs)