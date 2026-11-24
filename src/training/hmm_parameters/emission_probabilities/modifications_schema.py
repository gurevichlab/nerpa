from __future__ import annotations
from typing import (
    List,
    Tuple,
)
from enum import StrEnum
from src.general_type_aliases import Prob
import polars as pl

from src.antismash_parsing.bgc_variant_types import (
    BGC_Module,
    BGC_Module_Modification
)
from src.monomer_names_helper.monomer_names_helper import (
    Chirality,
    NRP_Monomer
)

class MODIFICATION_TYPE(StrEnum):
    METHYLATION = "METHYLATION"
    EPIMERIZATION = "EPIMERIZATION"

class ModificationsProbsSchema(StrEnum):
    MODIFICATION = "MODIFICATION"
    PRESENT_IN_BGC = "PRESENT_IN_BGC"
    PRESENT_IN_NRP = "PRESENT_IN_NRP"

    # P(PRESENT_IN_NRP | PRESENT_IN_BGC) for MODIFICATION
    PROBABILITY = "PROBABILITY"

    def _type(self) -> object:
        SCHEMA = ModificationsProbsSchema
        name_to_type = {
            SCHEMA.MODIFICATION: MODIFICATION_TYPE,
            SCHEMA.PRESENT_IN_BGC: bool,
            SCHEMA.PRESENT_IN_NRP: bool,
            SCHEMA.PROBABILITY: Prob,
        }

        assert name_to_type.keys() == {col.value for col in SCHEMA}, "Mismatch between enum members and type mapping"
        return name_to_type[self]


class ModificationsProbs_df:
    df: pl.DataFrame

    def __init__(self, df: pl.DataFrame):
        self.df = df

    @classmethod
    def from_match_emissions(
            cls,
            match_emissions: List[Tuple[BGC_Module, NRP_Monomer]]
    ) -> ModificationsProbs_df:
        MPS = ModificationsProbsSchema
        MT = MODIFICATION_TYPE

        rows: list[dict[str, object]] = []
        for bgc_module, nrp_monomer in match_emissions:
            bgc_mods = bgc_module.modifications
            bgc_meth = BGC_Module_Modification.METHYLATION in bgc_mods

            nrp_meth = nrp_monomer.methylated
            rows.append(
                {
                    MPS.MODIFICATION: MT.METHYLATION,
                    MPS.PRESENT_IN_BGC: bgc_meth,
                    MPS.PRESENT_IN_NRP: nrp_meth,                    
                }
            )

            bgc_chir = BGC_Module_Modification.EPIMERIZATION in bgc_mods
            if nrp_monomer.chirality == Chirality.UNKNOWN:
                continue
            else:
                nrp_chir = (
                    nrp_monomer.chirality == Chirality.D
                )

            rows.append(
                {
                    MPS.MODIFICATION: MT.EPIMERIZATION,
                    MPS.PRESENT_IN_BGC: bgc_chir,
                    MPS.PRESENT_IN_NRP: nrp_chir,
                }
            )

        df = pl.DataFrame(rows)
        df = df.group_by([
            MPS.MODIFICATION,
            MPS.PRESENT_IN_BGC,
            MPS.PRESENT_IN_NRP,
        ]).agg(
            pl.count().alias("COUNT"),
        )

        # P(PRESENT_IN_NRP | PRESENT_IN_BGC) per MODIFICATION
        df = (
            df
            .with_columns(
                pl.col("COUNT").sum()
                .over([MPS.MODIFICATION, MPS.PRESENT_IN_BGC])
                .alias("COUNT_TOTAL")
            )
            .with_columns(
                # Laplace pseudocounts
                ((pl.col("COUNT") + 1)
                 / (pl.col("COUNT_TOTAL") + 2))
                .cast(pl.Float64)
                .alias(MPS.PROBABILITY)
            )
            .drop("COUNT_TOTAL")
        )

        print("DEBUG INFO: ModificationsProbs_df.from_match_emissions()")
        print(df.head(5))

        return cls(df=df)
