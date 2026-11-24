from __future__ import annotations

import ast
from enum import StrEnum
from pathlib import Path
from typing import Any, Dict, List, Optional, NewType

import polars as pl


antiSMASH_MonomerName = NewType('antiSMASH_MonomerName', str)  # antiSMASH Short
NorineMonomerName = NewType('NorineMonomerName', str)
PARAS_MonomerName = NewType('PARAS_MonomerName', str)
NerpaResidue = NewType('NerpaResidue', str)

class MonomerNameFormat(StrEnum):
    ANTISMASH_SHORT = "antiSMASH_short"
    RBAN_NORINE = "rBAN/Norine"
    PARAS = "PARAS"


class MonomerType(StrEnum):
    NRPS = "NRPS"
    PKS = "PKS"
    PKS_HYBRID = "PKS_Hybrid"
    NOT_NRPS = "NOT_NRPS"
    NA = "NA"


class MonomerTag(StrEnum):
    ALPHA_HYDROXY_ACID = "alpha-hydroxy acid"


class MonomerNamesTableSchema(StrEnum):
    MONOMER_NAME = "MonomerName"
    NAME_FORMAT = "NameFormat"
    NERPA_RESIDUE = "NerpaResidue"
    MODIFICATIONS = "Modifications"
    TYPE = "Type"
    COMMENT = "Comment"
    TAGS = "Tags"

    def python_type(self) -> object:
        schema = MonomerNamesTableSchema
        column_types: Dict[MonomerNamesTableSchema, object] = {
            schema.MONOMER_NAME: str,
            schema.NAME_FORMAT: MonomerNameFormat,
            schema.NERPA_RESIDUE: NerpaResidue,
            schema.MODIFICATIONS: List[str],
            schema.TYPE: MonomerType,
            schema.COMMENT: Optional[str],
            schema.TAGS: Optional[List[MonomerTag]],
        }
        return column_types[self]

    def to_python_type(self, value: Any) -> Any:
        match self:
            case MonomerNamesTableSchema.MONOMER_NAME:
                return str(value)
            case MonomerNamesTableSchema.NAME_FORMAT:
                return MonomerNameFormat(value)
            case MonomerNamesTableSchema.NERPA_RESIDUE:
                return NerpaResidue(value)
            case MonomerNamesTableSchema.MODIFICATIONS:
                if isinstance(value, str):
                    return ast.literal_eval(value)
                return list(value)
            case MonomerNamesTableSchema.TYPE:
                return MonomerType(value)
            case MonomerNamesTableSchema.COMMENT:
                return str(value) if value is not None else None
            case MonomerNamesTableSchema.TAGS:
                tags = ast.literal_eval(value) if isinstance(value, str) else value
                return [MonomerTag(tag) for tag in tags]
            case _:
                raise ValueError(f"Unknown column: {self}")
            

    @classmethod
    def polars_schema(cls) -> Dict[str, pl.DataType]:
        column_types: Dict[MonomerNamesTableSchema, pl.DataType] = {
            cls.MONOMER_NAME: pl.Utf8,
            cls.NAME_FORMAT: pl.Utf8,
            cls.NERPA_RESIDUE: pl.Utf8,
            cls.MODIFICATIONS: pl.List(pl.Utf8),
            cls.TYPE: pl.Utf8,
            cls.COMMENT: pl.Utf8,
            cls.TAGS: pl.List(pl.Utf8),
        }
        return {
            column.value: data_type
            for column, data_type in column_types.items()
        }


class MonomerNamesTable_df:
    df: pl.DataFrame

    def __init__(self, df: pl.DataFrame):
        self.verify_schema(df)
        self.df = df

    @classmethod
    def from_tsv(cls, file_path: Path) -> MonomerNamesTable_df:
        schema = MonomerNamesTableSchema

        # Read serialized list columns as strings before parsing them.
        df = pl.read_csv(
            file_path,
            separator="\t",
            infer_schema=False,
        ).with_columns(
            pl.col(schema.MODIFICATIONS).map_elements(
                ast.literal_eval,
                return_dtype=pl.List(pl.Utf8),
            ),
            pl.col(schema.TAGS).map_elements(
                ast.literal_eval,
                return_dtype=pl.List(pl.Utf8),
            ),
        )

        return cls(df)

    @classmethod
    def verify_schema(cls, df: pl.DataFrame) -> None:
        expected_schema = MonomerNamesTableSchema.polars_schema()
        actual_schema = dict(df.schema)

        if actual_schema != expected_schema:
            raise ValueError(
                "DataFrame schema does not match the expected schema.\n"
                f"Expected: {expected_schema}\n"
                f"Actual: {actual_schema}"
            )

        for column in MonomerNamesTableSchema:
            for row_index, value in enumerate(df[column]):
                try:
                    column.to_python_type(value)
                except (TypeError, ValueError, SyntaxError) as error:
                    raise ValueError(
                        f"Invalid value in column {column.value!r}, "
                        f"row {row_index}: {value!r}"
                    ) from error
