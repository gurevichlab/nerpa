from __future__ import annotations
from typing import (
    Dict,
    Iterable,
    List,
    NamedTuple,
    Optional,
    NewType,
    Callable,
)
import pandas as pd
from enum import StrEnum

from src.antismash_parsing.antismash_parser_types import GeneId, BGC_ID
from src.antismash_parsing.bgc_variant_types import (
    BGC_Variant,
    BGC_Module,
    BGC_ID,
    BGC_Variant_ID,
)
    
from src.hmm.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    DetailedHMMStateType,
    StateIdx
)
from src.hmm.hmm_constructor.hmm_constructor_state_edge_context_relations import (
    filter_relevant_genomic_context
)
from dataclasses import dataclass, asdict

from src.monomer_names_helper import NerpaResidue
from src.rban_parsing.nrp_variant_types import (
    NRP_Variant,
    NRP_Variant_ID,
)
from src.rban_parsing.rban_monomer import rBAN_Monomer
from src.monomer_names_helper import NRP_Monomer
from src.antismash_parsing.genomic_context import ModuleGenomicContext
from src.testing.testing_types import CuratedAlignment
    
import polars as pl
from pathlib import Path
import numpy as np


class AlignmentWithBGCNRP(NamedTuple):
    curated_alignment: CuratedAlignment
    bgc_variant: BGC_Variant
    nrp_variant: NRP_Variant


def normalize_int_or_none(value: object) -> int | None:
    return (
        int(value)
        if value is not None
        else None
    )


class HMM_Chosen_Steps_Schema(StrEnum):
    BGC_ID = 'BGC_ID'  # BGC_Variant_ID
    NRP_ID = 'NRP_ID'  # NRP_Variant_ID
    STATE_TYPE = 'STATE_TYPE'  # DetailedHMMStateType
    EDGE_TYPE = 'EDGE_TYPE'  # DetailedHMMEdgeType
    MODULE_IDX = 'MODULE_IDX'  # int
    MODULE = 'MODULE'  # BGC_Module
    EMISSION = 'EMISSION'  # NRP_Monomer

    EDGE_FROM = 'EDGE_FROM'  # StateIdx
    EDGE_TO = 'EDGE_TO'  # StateIdx
    IDX_IN_INSERT_RUN = 'IDX_IN_INSERT_RUN'  # int

    # these columns are derived from MODULE and STATE_TYPE, and are not part of the original data
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
            SCHEMA.GENOMIC_CONTEXT: ModuleGenomicContext,
            SCHEMA.RELEVANT_GC: ModuleGenomicContext,
            SCHEMA.MODULE_IDX: int,
        }
        assert name_to_type.keys() == {col.value for col in SCHEMA}, "Mismatch between enum members and type mapping"
        return name_to_type[self]

    @classmethod
    def _polars_schema(
            cls,
            with_derived: bool = False
    ) -> Dict[str, pl.DataType]:
        SCHEMA = HMM_Chosen_Steps_Schema
        field_to_type = {
            SCHEMA.BGC_ID: pl.Object,
            SCHEMA.NRP_ID: pl.Object,
            SCHEMA.STATE_TYPE: pl.Utf8,  # store Enum.name
            SCHEMA.EDGE_TYPE: pl.Utf8,   # store Enum.name
            SCHEMA.MODULE_IDX: pl.Int64,
            SCHEMA.MODULE: pl.Object,
            SCHEMA.EMISSION: pl.Object,
            SCHEMA.EDGE_FROM: pl.Int64,
            SCHEMA.EDGE_TO: pl.Int64,
            SCHEMA.IDX_IN_INSERT_RUN: pl.Int64,
            SCHEMA.GENOMIC_CONTEXT: pl.Object,
            SCHEMA.RELEVANT_GC: pl.Object,
        }

        assert set(field_to_type.keys()) == set(SCHEMA), "Mismatch between enum members and type mapping"

        return {
            field.name: dtype
            for field, dtype in field_to_type.items()
            if with_derived or not field.is_derived()
        }

    def to_polars_dtype(self, value) -> object:
        SC = HMM_Chosen_Steps_Schema
        match self:
            case SC.STATE_TYPE | SC.EDGE_TYPE:
                return (
                    value.name
                    if value is not None
                    else None
                )
            case _:
                return value
        
        

    def is_derived(self) -> bool:
        # Columns not present in the original data;
        # computed from other columns
        SCHEMA = HMM_Chosen_Steps_Schema
        derived_columns = {
            SCHEMA.GENOMIC_CONTEXT,
            SCHEMA.RELEVANT_GC,
        }
        return self in derived_columns

    def tsv_str_converter(self) -> Callable[[object], str]:
        SC = HMM_Chosen_Steps_Schema
        def gc_representer(gc: ModuleGenomicContext) -> str:
            return str(sorted(f.name for f in gc))

        str_representer = {
            SC.BGC_ID: BGC_Variant_ID.to_str_short,
            SC.NRP_ID: NRP_Variant_ID.to_str_short,
            SC.STATE_TYPE: lambda x: x if x is not None else "",
            SC.EDGE_TYPE: lambda x: x if x is not None else "",
            SC.MODULE: lambda x: "",
            SC.EMISSION: (
                lambda mon: (
                    mon.to_base_mon().to_str_short()
                    if mon is not None
                    else ""
                )
            ),
            SC.EDGE_FROM: str,
            SC.EDGE_TO: str,
            SC.IDX_IN_INSERT_RUN: lambda x: str(x) if x is not None else "",
            SC.MODULE_IDX: lambda x: str(x) if x is not None else "",
            SC.GENOMIC_CONTEXT: (
                lambda x: (
                    gc_representer(x)
                    if x is not None
                    else "")
            ),
            SC.RELEVANT_GC: (
                lambda x: (
                    gc_representer(x)
                    if x is not None
                    else "")
            ),
        }
        return str_representer[self]

    @classmethod
    def missing_columns_in_df(cls, df: pl.DataFrame) -> List[HMM_Chosen_Steps_Schema]:
        return [
            col
            for col in HMM_Chosen_Steps_Schema
            if (col not in df.columns
                and not col.is_derived())
        ]



class HMM_Chosen_Steps_df:
    df: pl.DataFrame

    def _fill_derived_columns(self):
        SC = HMM_Chosen_Steps_Schema

        genomic_context: List[ModuleGenomicContext | None] = [
            m.genomic_context if m is not None else None
            for m in self.df[SC.MODULE]
        ]

        relevant_gc: List[ModuleGenomicContext | None] = [
            (
                filter_relevant_genomic_context(
                    gc,
                    DetailedHMMStateType[state_type]
                )
                if gc is not None
                else None
            )
            for state_type, gc in zip(self.df[SC.STATE_TYPE], genomic_context)
        ]

        self.df = self.df.with_columns(
            pl.Series(SC.GENOMIC_CONTEXT, genomic_context, dtype=pl.Object),
            pl.Series(SC.RELEVANT_GC, relevant_gc, dtype=pl.Object),
        )

    def __init__(
            self,
            df: pl.DataFrame,
            derive_columns: bool = True
    ):
        if (missing_cols := HMM_Chosen_Steps_Schema.missing_columns_in_df(df)):
            raise ValueError(
                f"HMM_Chosen_Steps_df DataFrame is missing required columns: "
                f"{missing_cols}"
            )

        self.df = df
        if derive_columns:
            self._fill_derived_columns()

    @classmethod
    def from_rows(
            cls,
            rows: List[Dict[HMM_Chosen_Steps_Schema, object]],
            derive_columns: bool = True
    ) -> HMM_Chosen_Steps_df:
        # Polars' row-wise constructor ("from dicts") will normalize dataclasses /
        # namedtuple-ish objects into plain Python dicts. Build column-wise to
        # preserve exact Python objects in pl.Object columns.
        from enum import Enum, IntEnum, StrEnum
        SC = HMM_Chosen_Steps_Schema
        schema: Dict[str, pl.DataType] = (
            SC._polars_schema(with_derived=not derive_columns)
        )
        column_data: Dict[str, List[object]] = {
            col.name: [
                col.to_polars_dtype(row[col])
                for row in rows
            ]
            for col in SC
        }

        print(
            "Num gc context == None: "
            f"{len([v for v in column_data[SC.GENOMIC_CONTEXT.name] if v is None])}"
        )

        print(
            "Checking schema match\n"
            f"column_data-schema:\n{column_data.keys() - schema.keys()}\n\n"
            f"schema-column_data:\n{schema.keys() - column_data.keys()}\n"
            f"Keys match: {schema.keys() == column_data.keys()}"
        )

        df = pl.DataFrame(
            column_data,
            schema=schema
        )

        return HMM_Chosen_Steps_df(df, derive_columns=derive_columns)

    def to_tsv_rows(self) -> List[Dict[str, str]]:
        rows: List[Dict[str, str]] = []
        for row in self.df.iter_rows(named=True):
            row_dict: Dict[str, str] = {}
            for col in HMM_Chosen_Steps_Schema:
                value = row[col]
                str_representer = col.tsv_str_converter()
                row_dict[col] = str_representer(value)
            rows.append(row_dict)
        return rows

    def write_tsv(self, file_path: Path) -> None:
        rows = self.to_tsv_rows()
        df = pl.DataFrame(rows)
        df.write_csv(file_path, separator='\t')

    def transition_cnts(self) -> pl.DataFrame:
        SC = HMM_Chosen_Steps_Schema
        df = (
            self.df
            .select([SC.STATE_TYPE, SC.EDGE_TYPE, SC.RELEVANT_GC])
            .group_by([SC.STATE_TYPE, SC.EDGE_TYPE, SC.RELEVANT_GC])
            .agg(pl.count().alias('count'))
        )
        return df

    def write_transition_cnts_tsv(self, file_path: Path) -> None:
        SC = HMM_Chosen_Steps_Schema
        rows = self.to_tsv_rows()
        df = (
            pl.DataFrame(rows)
            .select([SC.STATE_TYPE, SC.EDGE_TYPE, SC.RELEVANT_GC])
            .group_by([SC.STATE_TYPE, SC.EDGE_TYPE, SC.RELEVANT_GC])
            .agg(pl.count().alias('count'))
            .sort([SC.STATE_TYPE, SC.EDGE_TYPE, SC.RELEVANT_GC])
        )
        
        df.write_csv(file_path, separator='\t')
        


class SuccessCnts(NamedTuple):
    failures: int
    successes: int

    def __add__(self, other: SuccessCnts) -> SuccessCnts:
        if not isinstance(other, SuccessCnts):
            return NotImplemented
        return SuccessCnts(
            failures=self.failures + other.failures,
            successes=self.successes + other.successes
        )

    def total(self) -> int:
        return self.failures + self.successes

    def estimate_prob(self, pseudocounts: bool = False) -> Prob:
        if pseudocounts:
            return (self.successes + 1) / (self.total() + 2)
        if self.total() == 0:
            raise ValueError('Cannot estimate probability with zero total counts.')
        return self.successes / self.total()
