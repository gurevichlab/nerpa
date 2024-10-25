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
from auxilary_types import StepInfo




def get_matches_stats(match: Match,
                      bgc_variant: BGC_Variant,
                      nrp_variant: NRP_Variant) -> MatchStats:
    module_skips, gene_skips, fragment_skips = get_skips(match, bgc_variant)
    return MatchStats(match_steps=match.alignment)

def get_matches_stats_for_bgc(matches_with_bgc_nrp: Iterable[MatchWithBGCNRP]):
    pass