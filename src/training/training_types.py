from typing import (
    Dict,
    List,
    NamedTuple,
    Optional
)

from src.matching.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    GenomicContext,
    EdgeKey
)
from src.matching.match_type import Match
from src.data_types import (
    BGC_Module,
    BGC_Variant,
    NRP_Monomer,
    NRP_Variant
)
from dataclasses import dataclass


class MatchWithBGCNRP(NamedTuple):
    match: Match
    bgc_variant: BGC_Variant
    nrp_variant: Optional[NRP_Variant]


class ChoicesCnts(NamedTuple):
    NOT_CHOSEN: int
    CHOSEN: int


class MatchEmissionInfo(NamedTuple):
    bgc_module: BGC_Module
    nrp_monomer: NRP_Monomer


@dataclass
class DataForTraining:
    edge_choices_cnts: Dict[DetailedHMMEdgeType, Dict[Optional[GenomicContext], ChoicesCnts]]
    match_emissions: List[MatchEmissionInfo]  # (module, monomer)
    # insert_emissions: List[Tuple[BGC_Module, NRP_Monomer]]
    # insert_at_start_emissions: List[Tuple[BGC_Module, NRP_Monomer]]


# for the case of several insertions in a row I want to distinguish them
class ExtendedEdgeKey(NamedTuple):
    edge_key: EdgeKey
    num_insertions: int = 0


class EdgeInfo(NamedTuple):
    edge_type: DetailedHMMEdgeType
    genomic_context: GenomicContext


@dataclass
class PathTurnInfo:
    chosen_edge_key: ExtendedEdgeKey
    chosen_edge_info: EdgeInfo
    other_edges_info: List[EdgeInfo]
    emission_info: Optional[MatchEmissionInfo] = None
