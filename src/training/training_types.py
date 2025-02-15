from typing import NamedTuple, Optional
from src.matching.matching_types_match import Match
from src.data_types import BGC_Variant, NRP_Variant


class MatchWithBGCNRP(NamedTuple):
    match: Match
    bgc_variant: BGC_Variant
    nrp_variant: Optional[NRP_Variant]


