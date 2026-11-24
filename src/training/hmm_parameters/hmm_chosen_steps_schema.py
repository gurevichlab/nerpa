from __future__ import annotations
from typing import (
    Callable,
    Dict,
    List,
    NamedTuple,
    Optional, NewType
)
import pandas as pd
from enum import Enum

from src.antismash_parsing.antismash_parser_types import GeneId, BGC_ID
from src.antismash_parsing.bgc_variant_types import (
    BGC_Variant,
    BGC_Module,
    BGC_ID,
    BGC_Variant_ID,
)
    
from src.hmm.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    GenomicContext,
    DetailedHMMStateType, StateIdx
)
from src.hmm.hmm_constructor.hmm_constructor_state_edge_context_relations import (
    filter_relevant_genomic_context
)
from dataclasses import dataclass, asdict

from src.monomer_names_helper.monomer_names_helper import NerpaResidue
from src.rban_parsing.nrp_variant_types import (
    NRP_Variant,
    NRP_Variant_ID,
)
from src.rban_parsing.rban_monomer import rBAN_Monomer
from src.monomer_names_helper.monomer_names_helper import NRP_Monomer
from src.testing.testing_types import CuratedAlignment
import polars as pl


class HMM_Chosen_Steps_Schema(Enum):
    BGC_ID = 'BGC_ID'  # BGC_Variant_ID
    NRP_ID = 'NRP_ID'  # NRP_Variant_ID
    STATE_TYPE = 'STATE_TYPE'  # DetailedHMMStateType
    EDGE_TYPE = 'EDGE_TYPE'  # DetailedHMMEdgeType
    MODULE = 'MODULE'  # BGC_Module
    EMISSION = 'EMISSION'  # NRP_Monomer

    EDGE_FROM = 'EDGE_FROM'  # StateIdx
    EDGE_TO = 'EDGE_TO'  # StateIdx
    IDX_IN_INSERT_RUN = 'IDX_IN_INSERT_RUN'  # int

    # these columns are derived from MODULE and STATE_TYPE, and are not part of the original data
    MODULE_IDX = 'MODULE_IDX'  # int
    GENOMIC_CONTEXT = 'GENOMIC_CONTEXT'  # GenomicContext
    RELEVANT_GC = 'RELEVANT_GC'  # GenomicContext

    def _type(self) -> object:
        SCHEMA = HMM_Chosen_Steps_Schema
        name_to_type = {
            SCHEMA.BGC_ID: BGC_Variant_ID,
            SCHEMA.NRP_ID: NRP_Variant_ID,
            SCHEMA.STATE_TYPE: DetailedHMMStateType,
            SCHEMA.EDGE_TYPE: DetailedHMMEdgeType,
            SCHEMA.MODULE: BGC_Module,
            SCHEMA.EMISSION: NRP_Monomer,
            SCHEMA.EDGE_FROM: StateIdx,
            SCHEMA.EDGE_TO: StateIdx,
            SCHEMA.IDX_IN_INSERT_RUN: int,

            SCHEMA.GENOMIC_CONTEXT: GenomicContext,
            SCHEMA.RELEVANT_GC: GenomicContext,
            SCHEMA.MODULE_IDX: int,
        }

        assert name_to_type.keys() == {col.value for col in SCHEMA}, "Mismatch between enum members and type mapping"
        return name_to_type[self]

    def is_derived(self) -> bool:
        # Columns not present in the original data;
        # computed from other columns
        SCHEMA = HMM_Chosen_Steps_Schema
        derived_columns = {
            SCHEMA.MODULE_IDX,
            SCHEMA.GENOMIC_CONTEXT,
            SCHEMA.RELEVANT_GC,
        }
        return self in derived_columns

        

    def tsv_str_converter(self) -> Callable:
        SC = HMM_Chosen_Steps_Schema
        def gc_representer(gc: ModuleGenomicContext) -> str:
            return str(gc)

        str_representer = {
            SC.BGC_ID: BGC_Variant_ID.to_str_short,
            SC.NRP_ID: NRP_Variant_ID.to_str_short,
            SC.STATE_TYPE: lambda x: x.name,
            SC.EDGE_TYPE: lambda x: x.name,
            SC.MODULE: lambda x: "",
            SC.EMISSION: (
                lambda mon: (
                    mon.to_str_short()
                    if mon is not None
                    else ""
                )
            ),
            SC.EDGE_FROM: str,
            SC.EDGE_TO: str,
            SC.IDX_IN_INSERT_RUN: lambda x: str(x) if x is not None else "",
            SC.MODULE_IDX: lambda x: str(x) if x is not None else "",
            SC.GENOMIC_CONTEXT: lambda x: str(x) if x is not None else "",
            SC.RELEVANT_GC: lambda x: str(x) if x is not None else "",
        }
        raise NotImplementedError("TSV string converter is not implemented yet.")
