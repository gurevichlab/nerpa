from typing import NamedTuple
from src.matching.matching_types_match import Match
from src.data_types import BGC_Variant, NRP_Variant
from src.matching.matching_types_alignment_step import AlignmentStep
from abc import ABC, abstractmethod


class MatchWithBGCNRP(NamedTuple):
    match: Match
    bgc_variant: BGC_Variant
    nrp_variant: NRP_Variant


