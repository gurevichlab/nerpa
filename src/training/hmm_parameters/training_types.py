from __future__ import annotations
from typing import (
    Dict,
    List,
    NamedTuple,
    Optional
)

from src.antismash_parsing.antismash_parser_types import GeneId
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
    NRP_Variant, BGC_ID
)
from dataclasses import dataclass

from src.monomer_names_helper import MonomerResidue
from src.rban_parsing.rban_monomer import rBAN_Monomer


class MatchWithBGCNRP(NamedTuple):
    match: Match
    bgc_variant: BGC_Variant
    nrp_variant: Optional[NRP_Variant]


class ChoicesCnts(NamedTuple):
    NOT_CHOSEN: int
    CHOSEN: int

    def to_prob(self, pseudocounts: bool = False) -> float:
        if pseudocounts:
            return (self.CHOSEN + 1) / (self.NOT_CHOSEN + self.CHOSEN + 2)
        else:
            return self.CHOSEN / (self.NOT_CHOSEN + self.CHOSEN)


class MatchEmissionInfo(NamedTuple):
    bgc_id: BGC_ID
    bgc_module: BGC_Module
    nrp_monomer: rBAN_Monomer


class MatchEmissionKey(NamedTuple):
    bgc_id: BGC_ID
    gene_id: GeneId
    a_domain_idx: int
    resisue: MonomerResidue

    @classmethod
    def from_match_emission_info(cls, match_emission_info: MatchEmissionInfo) -> MatchEmissionKey:
        return cls(bgc_id=match_emission_info.bgc_id,
                   gene_id=match_emission_info.bgc_module.gene_id,
                   a_domain_idx=match_emission_info.bgc_module.a_domain_idx,
                   resisue=match_emission_info.nrp_monomer.residue)


@dataclass
class DataForTraining:
    edge_choices_cnts: Dict[DetailedHMMEdgeType, Dict[Optional[GenomicContext], ChoicesCnts]]
    match_emissions: List[MatchEmissionInfo]  # (module, monomer)
    # insert_emissions: List[Tuple[BGC_Module, NRP_Monomer]]
    # insert_at_start_emissions: List[Tuple[BGC_Module, NRP_Monomer]]


# for the case of several insertions in a row I want to distinguish them
class ExtendedEdgeKey(NamedTuple):
    edge_type: DetailedHMMEdgeType
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
