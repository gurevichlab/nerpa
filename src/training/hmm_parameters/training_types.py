from __future__ import annotations
from typing import (
    Dict,
    List,
    NamedTuple,
    Optional, NewType
)
import pandas as pd

from src.antismash_parsing.antismash_parser_types import GeneId, BGC_ID
from src.antismash_parsing.bgc_variant_types import BGC_Variant, BGC_Module
from src.hmm.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    GenomicContext,
    DetailedHMMStateType, StateIdx
)
from src.hmm.hmm_constructor.hmm_constructor_state_edge_context_relations import (
    filter_relevant_genomic_context
)
from dataclasses import dataclass, asdict

from src.monomer_names_helper import NerpaResidue
from src.rban_parsing.nrp_variant_types import NRP_Variant
from src.rban_parsing.rban_monomer import rBAN_Monomer
from src.testing.testing_types import CuratedAlignment
import polars as pl


class AlignmentWithBGCNRP(NamedTuple):
    curated_alignment: CuratedAlignment
    bgc_variant: BGC_Variant
    nrp_variant: NRP_Variant


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
    resisue: NerpaResidue
    state_type: DetailedHMMStateType

    @classmethod
    def from_emission_info(cls, match_emission_info: EmissionInfo) -> EmissionKey:
        return cls(bgc_id=match_emission_info.bgc_id,
                   gene_id=match_emission_info.bgc_module.gene_id if match_emission_info.bgc_module else GeneId(''),
                   a_domain_idx=match_emission_info.bgc_module.a_domain_idx if match_emission_info.bgc_module else -1,
                   resisue=match_emission_info.nrp_monomer.residue,
                   state_type=match_emission_info.state_type)


class HMM_Chosen_Steps_Schema:
    BGC_ID = 'BGC_ID'  # BGC_Variant_ID
    NRP_ID = 'NRP_ID'  # NRP_Variant_ID
    STATE_TYPE = 'STATE_TYPE'  # DetailedHMMStateType
    EDGE_TYPE = 'EDGE_TYPE'  # DetailedHMMEdgeType
    MODULE = 'MODULE'  # BGC_Module
    EMISSION = 'EMISSION'  # NRP_Monomer

    EDGE_FROM = 'EDGE_FROM'  # StateIdx
    EDGE_TO = 'EDGE_TO'  # StateIdx
    NUM_INSERTIONS = 'NUM_INSERTIONS'  # int

    # these columns are derived from MODULE and STATE_TYPE, and are not part of the original data
    GENOMIC_CONTEXT = 'GENOMIC_CONTEXT'  # GenomicContext
    RELEVANT_GC = 'RELEVANT_GC'  # GenomicContext


class HMM_Chosen_Steps_df(pl.DataFrame):
    @classmethod
    def from_rows(cls, rows: List[Dict[str, object]]) -> HMM_Chosen_Steps_df:
        ECS = HMM_Chosen_Steps_Schema

        df = pl.DataFrame(rows)
        missing_cols = [
            col
            for col in HMM_Chosen_Steps_Schema.__dict__.values()
            if col not in df.columns
        ]

        if missing_cols:
            raise ValueError(f'HMM_Chosen_Steps_df DataFrame is missing required columns: {missing_cols}')


        df = df.with_columns(
            pl.col(ECS.MODULE)
            .map_elements(
                lambda m: m.genomic_context if m is not None else None,
                return_dtype=pl.Object,
            )
            .alias(ECS.GENOMIC_CONTEXT)
        )

        df = df.with_columns(
            pl.struct([ECS.GENOMIC_CONTEXT, ECS.STATE_TYPE])
            .map_elements(
                lambda s: filter_relevant_genomic_context(s[ECS.STATE_TYPE], s[ECS.GENOMIC_CONTEXT]),
                return_dtype=pl.Object,
            )
            .alias(ECS.RELEVANT_GC)
        )

        return HMM_Chosen_Steps_df(df)


@dataclass
class DataForTraining:
    edge_choices_per_state_cntxt: Dict[DetailedHMMStateType, Dict[GenomicContext, Dict[DetailedHMMEdgeType, int]]]
    edge_choices_wo_filtering: Dict[DetailedHMMStateType, Dict[GenomicContext, Dict[DetailedHMMEdgeType, int]]]
    edge_choices_df: HMM_Chosen_Steps_df
    #edge_choices_cnts: Dict[DetailedHMMEdgeType, Dict[GenomicContext, ChoicesCnts]]
    emissions: List[EmissionInfo]  # (module, monomer)
    #chosen_edges_occurrences: Dict[DetailedHMMEdgeType, Dict[GenomicContext, List[BGC_ID]]]
    # insert_emissions: List[Tuple[BGC_Module, NRP_Monomer]]
    # insert_at_start_emissions: List[Tuple[BGC_Module, NRP_Monomer]]

    # def to_dict(self) -> Dict[str, object]:
    #     def context_to_str(context: GenomicContext) -> str:
    #         return ('('
    #                 + ','.join(context_feature.name
    #                            for context_feature in context)
    #                 + ')')
    #     return {
    #         'edge_choices_cnts': {edge_type.name: {context_to_str(context): cnts._asdict()
    #                                                for context, cnts in context_cnts.items()}
    #                               for edge_type, context_cnts in self.edge_choices_cnts.items()},
    #         'emissions': [emission.to_dict() for emission in self.emissions],
    #         'chosen_edges_occurrences': {edge_type.name: {context_to_str(context): bgc_ids
    #                                                       for context, bgc_ids in context_bgc_ids.items()}
    #                                      for edge_type, context_bgc_ids in self.chosen_edges_occurrences.items()}
    #     }

    def to_dict(self) -> Dict[str, object]:
        def context_to_str(context: GenomicContext) -> str:
            feature_names = [context_feature.name for context_feature in context]
            if len(feature_names) == 1:
                return f'({feature_names[0]},)'
            else:
                return f'({",".join(feature_names)})'

        return {
            'edge_choices_per_state_ctxt': {
                state_type.name: {
                    context_to_str(context): dict(cnts)
                    for context, cnts in context_cnts.items()
                }
                for state_type, context_cnts in self.edge_choices_per_state_cntxt.items()
            },
            'emissions': [emission.to_dict() for emission in self.emissions],
        }


class EdgeKey(NamedTuple):
    start: StateIdx
    end: StateIdx
    num_insertions: int = 0

class EdgeInfo(NamedTuple):
    edge_type: DetailedHMMEdgeType
    edge_start_type: DetailedHMMStateType
    genomic_context: GenomicContext


@dataclass
class PathTurnInfo:
    bgc_id: BGC_ID
    chosen_edge_key: EdgeKey
    chosen_edge_info: EdgeInfo
    other_edges_info: List[EdgeInfo]
