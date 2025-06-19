from __future__ import annotations
from typing import (
    Dict,
    List,
    NamedTuple,
    Optional, Literal
)

from src.antismash_parsing.antismash_parser_types import GeneId
from src.matching.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    GenomicContext,
    EdgeKey, DetailedHMMStateType
)
from src.matching.match_type import Match
from src.data_types import (
    BGC_Module,
    BGC_Variant,
    NRP_Monomer,
    NRP_Variant, BGC_ID
)
from dataclasses import dataclass, asdict

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


class EmissionInfo(NamedTuple):
    bgc_id: BGC_ID
    bgc_module: Optional[BGC_Module]
    nrp_monomer: rBAN_Monomer
    state_type: DetailedHMMStateType

    def to_dict(self) -> Dict[str, object]:
        return {
            'bgc_id': self.bgc_id._asdict(),
            'bgc_module': asdict(self.bgc_module) if self.bgc_module is not None else None,
            'nrp_monomer': self.nrp_monomer._asdict(),
            'state_type': self.state_type.name
        }


class EmissionKey(NamedTuple):
    bgc_id: BGC_ID
    gene_id: GeneId
    a_domain_idx: int
    resisue: MonomerResidue
    state_type: DetailedHMMStateType

    @classmethod
    def from_emission_info(cls, match_emission_info: EmissionInfo) -> EmissionKey:
        return cls(bgc_id=match_emission_info.bgc_id,
                   gene_id=match_emission_info.bgc_module.gene_id if match_emission_info.bgc_module else GeneId(''),
                   a_domain_idx=match_emission_info.bgc_module.a_domain_idx if match_emission_info.bgc_module else -1,
                   resisue=match_emission_info.nrp_monomer.residue,
                   state_type=match_emission_info.state_type)


@dataclass
class DataForTraining:
    edge_choices_cnts: Dict[DetailedHMMEdgeType, Dict[GenomicContext, ChoicesCnts]]
    emissions: List[EmissionInfo]  # (module, monomer)
    chosen_edges_occurrences: Dict[DetailedHMMEdgeType, Dict[GenomicContext, List[BGC_ID]]]
    # insert_emissions: List[Tuple[BGC_Module, NRP_Monomer]]
    # insert_at_start_emissions: List[Tuple[BGC_Module, NRP_Monomer]]

    def to_dict(self) -> Dict[str, object]:
        def context_to_str(context: GenomicContext) -> str:
            return ('('
                    + ','.join(context_feature.name
                               for context_feature in context)
                    + ')')
        return {
            'edge_choices_cnts': {edge_type.name: {context_to_str(context): cnts._asdict()
                                                   for context, cnts in context_cnts.items()}
                                  for edge_type, context_cnts in self.edge_choices_cnts.items()},
            'emissions': [emission.to_dict() for emission in self.emissions],
            'chosen_edges_occurrences': {edge_type.name: {context_to_str(context): bgc_ids
                                                          for context, bgc_ids in context_bgc_ids.items()}
                                         for edge_type, context_bgc_ids in self.chosen_edges_occurrences.items()}
        }


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
    bgc_id: BGC_ID
    chosen_edge_key: ExtendedEdgeKey
    chosen_edge_info: EdgeInfo
    other_edges_info: List[EdgeInfo]
