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
from src.general_type_aliases import Prob, LogProb
from dataclasses import dataclass, asdict

from src.monomer_names_helper import NerpaResidue
from src.rban_parsing.nrp_variant_types import (
    NRP_Variant,
    NRP_Variant_ID,
)
from src.rban_parsing.rban_monomer import rBAN_Monomer
from src.monomer_names_helper import NRP_Monomer
from src.antismash_parsing.genomic_context import (
    ModuleGenomicContext,
    ModuleGenomicContextFeature
)
from src.testing.testing_types import CuratedAlignment
from src.generic.numeric import safe_log
    
import polars as pl
from pathlib import Path
import numpy as np

# can't import short_aliases here because of circular import, so we define them here
MCF = ModuleGenomicContextFeature
GC = ModuleGenomicContext
ST = DetailedHMMStateType
ET = DetailedHMMEdgeType


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
    AA34 = 'AA34'  # str

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
            SCHEMA.GENOMIC_CONTEXT: pl.List(pl.Utf8),
            SCHEMA.RELEVANT_GC: pl.List(pl.Utf8),
            SCHEMA.AA34: pl.Utf8,
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
            case SC.GENOMIC_CONTEXT | SC.RELEVANT_GC:
                return sorted(f.name for f in value)
            case _:
                return value

    def from_polars_dtype(self, value) -> object:
        SC = HMM_Chosen_Steps_Schema
        match self:
            case SC.STATE_TYPE:
                return (
                    DetailedHMMStateType[value]
                    if value is not None
                    else None
                )
            case SC.EDGE_TYPE:
                return (
                    DetailedHMMEdgeType[value]
                    if value is not None
                    else None
                )
            case SC.GENOMIC_CONTEXT | SC.RELEVANT_GC:
                return (
                    frozenset({ModuleGenomicContextFeature[f] for f in value})
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
            SCHEMA.AA34,
        }
        return self in derived_columns

    def tsv_str_converter(self) -> Callable[[object], str]:
        SC = HMM_Chosen_Steps_Schema

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
                lambda gc: (
                    '[' + ', '.join(gc) + ']'
                    if gc is not None
                    else ""
                )
            ),
            SC.RELEVANT_GC: (
                lambda gc: (
                    '[' + ', '.join(gc) + ']'
                    if gc is not None
                    else ""
                )
            ),
            SC.AA34: lambda x: x if x is not None else "",
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

        genomic_context: List[ModuleGenomicContext] = [
            (m.genomic_context if m is not None else frozenset())
            for m in self.df[SC.MODULE]
        ]

        relevant_gc: List[ModuleGenomicContext] = [
            filter_relevant_genomic_context(
                gc,
                DetailedHMMStateType[state_type]
            )
            for state_type, gc in zip(self.df[SC.STATE_TYPE], genomic_context)
        ]

        genomic_context_plain: List[List[str]] = [
            [f.name for f in gc]
            for gc in genomic_context
        ]
        relevant_gc_plain: List[List[str]] = [
            [f.name for f in gc]
            for gc in relevant_gc
        ]

        aa34s: List[str | None] = [
            m.aa34_code if m is not None else None
            for m in self.df[SC.MODULE]
        ]

        self.df = self.df.with_columns(
            pl.Series(
                SC.GENOMIC_CONTEXT,
                genomic_context_plain,
                dtype=pl.List(pl.Utf8)
            ),
            pl.Series(
                SC.RELEVANT_GC,
                relevant_gc_plain,
                dtype=pl.List(pl.Utf8)
            ),
            pl.Series(
                SC.AA34,
                aa34s,
                dtype=pl.Utf8
            )
        )

    def __init__(
            self,
            df: pl.DataFrame | None = None,
            derive_columns: bool = True
    ):
        if df is None:
            # Create an empty DataFrame with the expected schema
            SC = HMM_Chosen_Steps_Schema
            df = pl.DataFrame(
                {
                    col.name: (
                        pl.Series(
                            col.name,
                            [],
                            dtype=SC._polars_schema()[col.name]
                        )
                    )
                    for col in SC
                }
            )
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
        if not rows:
            return HMM_Chosen_Steps_df()

        SC = HMM_Chosen_Steps_Schema
        schema: Dict[str, pl.DataType] = (
            SC._polars_schema(with_derived=not derive_columns)
        )
        column_data: Dict[str, List[object]] = {
            col.name: [
                col.to_polars_dtype(row[col])
                for row in rows
            ]
            for col in rows[0].keys()
        }


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


class Estimated_Transition_Probs_Schema(StrEnum):
    STATE_TYPE = 'STATE_TYPE'  # DetailedHMMStateType
    RELEVANT_GC = 'RELEVANT_GC'  # GenomicContext
    EDGE_TYPE = 'EDGE_TYPE'  # DetailedHMMEdgeType
    PROB = 'PROB'  # Prob
    COMMENTS = 'COMMENTS'  # List[str]

    def to_polars_dtype(
            self,
            value: object,
    ) -> object:
        match self:
            case Estimated_Transition_Probs_Schema.STATE_TYPE:
                return value.name if value is not None else None
            case Estimated_Transition_Probs_Schema.RELEVANT_GC:
                return sorted(f.name for f in value) if value is not None else None
            case Estimated_Transition_Probs_Schema.EDGE_TYPE:
                return value.name if value is not None else None
            case Estimated_Transition_Probs_Schema.PROB:
                return float(value) if value is not None else None
            case Estimated_Transition_Probs_Schema.COMMENTS:
                return sorted(value) if value is not None else None
            case _:
                raise ValueError(f"Unknown column: {col}")
    
class ProbabilityEstimate(NamedTuple):
    prob: Prob
    comments: list[str]

class Estimated_Transition_Probs:
    df: pl.DataFrame

    def __init__(self, df: pl.DataFrame | None = None):
        if df is not None:
            self.df = df
            return

        # Ensure an empty instance still has the expected columns + dtypes.
        SC = Estimated_Transition_Probs_Schema
        self.df = pl.DataFrame(
            {
                SC.STATE_TYPE: pl.Series(SC.STATE_TYPE, [], dtype=pl.Utf8),
                SC.RELEVANT_GC: pl.Series(SC.RELEVANT_GC, [], dtype=pl.List(pl.Utf8)),
                SC.EDGE_TYPE: pl.Series(SC.EDGE_TYPE, [], dtype=pl.Utf8),
                SC.PROB: pl.Series(SC.PROB, [], dtype=pl.Float64),
                SC.COMMENTS: pl.Series(SC.COMMENTS, [], dtype=pl.List(pl.Utf8)),
            }
        )

    @classmethod
    def from_rows(
            cls,
            rows: List[Dict[Estimated_Transition_Probs_Schema, object]]
    ) -> Estimated_Transition_Probs:
        SC = Estimated_Transition_Probs_Schema
        schema: Dict[str, pl.DataType] = {
            SC.STATE_TYPE: pl.Utf8,
            SC.RELEVANT_GC: pl.List(pl.Utf8),
            SC.EDGE_TYPE: pl.Utf8,
            SC.PROB: pl.Float64,
            SC.COMMENTS: pl.List(pl.Utf8),
        }
        column_data: Dict[str, List[object]] = {
            col.name: [
                col.to_polars_dtype(row[col])
                for row in rows
            ]
            for col in SC
        }

        df = pl.DataFrame(
            column_data,
            schema=schema
        )

        return Estimated_Transition_Probs(df)


    def to_dict(
            self,
            log_probs: bool = False
    ) -> Dict[ST, Dict[GC, Dict[ET, Prob | LogProb]]]:
        SC = Estimated_Transition_Probs_Schema
        transition_probs: Dict[ST, Dict[GC, Dict[ET, Prob]]] = {}
        for row in self.df.iter_rows(named=True):
            state_type: ST = ST[row[SC.STATE_TYPE]]
            relevant_gc: GC = frozenset({ModuleGenomicContextFeature[f] for f in row[SC.RELEVANT_GC]})
            edge_type: ET = ET[row[SC.EDGE_TYPE]]
            prob: Prob = row[SC.PROB]
            if log_probs:
                prob = safe_log(prob)

            if state_type not in transition_probs:
                transition_probs[state_type] = {}
            if relevant_gc not in transition_probs[state_type]:
                transition_probs[state_type][relevant_gc] = {}
            transition_probs[state_type][relevant_gc][edge_type] = prob

        return transition_probs

    @classmethod
    def from_dict(
            cls,
            transition_probs: Dict[ST, Dict[GC, Dict[ET, Prob | ProbabilityEstimate]]]
    ) -> Estimated_Transition_Probs:
        ETPS = Estimated_Transition_Probs_Schema
        rows: List[Dict[ETPS, object]] = []
        for state_type, gc_dict in transition_probs.items():
            for relevant_gc, edge_dict in gc_dict.items():
                for edge_type, prob_est in edge_dict.items():
                    if isinstance(prob_est, ProbabilityEstimate):
                        prob = prob_est.prob
                        comments = prob_est.comments
                    else:
                        prob = prob_est
                        comments = []
                    rows.append({
                        ETPS.STATE_TYPE: state_type,
                        ETPS.RELEVANT_GC: relevant_gc,
                        ETPS.EDGE_TYPE: edge_type,
                        ETPS.PROB: prob,
                        ETPS.COMMENTS: comments,
                    })
        return cls.from_rows(rows)
                    
    def update(self, other: Estimated_Transition_Probs) -> None:
        self.df = pl.concat([self.df, other.df], how='vertical').unique()

    def write_tsv(self, file_path: Path) -> None:
        SC = Estimated_Transition_Probs_Schema

        # Polars writes List columns in a non-TSV-friendly way; stringify them.
        df_to_write = self.df.with_columns(
            pl.when(pl.col(SC.RELEVANT_GC).is_null())
            .then(pl.lit(""))
            .otherwise(
                pl.concat_str(
                    [
                        pl.lit("["),
                        pl.col(SC.RELEVANT_GC).list.join(", "),
                        pl.lit("]"),
                    ]
                )
            )
            .alias(SC.RELEVANT_GC),
            pl.when(pl.col(SC.COMMENTS).is_null())
            .then(pl.lit(""))
            .otherwise(
                pl.concat_str(
                    [
                        pl.lit("['"),
                        pl.col(SC.COMMENTS).list.join("', '"),
                        pl.lit("']"),
                    ]
                )
            )
            .alias(SC.COMMENTS),
        )

        df_to_write.write_csv(file_path, separator="\t")

