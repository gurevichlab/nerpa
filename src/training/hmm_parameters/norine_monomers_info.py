from __future__ import annotations
from typing import (
    Dict,
    List,
    Tuple,
)
from dataclasses import dataclass
from src.monomer_names_helper.monomer_names_helper import NerpaResidue
import yaml
from pathlib import Path
from enum import StrEnum
import polars as pl
from src.general_type_aliases import Prob
from src.training.hmm_parameters.emission_probabilities.modifications_schema import MODIFICATION_TYPE
from src.monomer_names_helper.monomer_names_helper import (
    Chirality,
    UNKNOWN_RESIDUE
)
from collections import Counter


class NorineMonomersInfoSchema(StrEnum):
    NORINE_NAME = "Norine_name"
    COUNT = "count"
    NERPA_RESIDUE = "Nerpa_residue"
    METHYLATED = "methylated"
    CHIRALITY = "chirality"
    IS_PKS_HYBRID = "is_pks_hybrid"


def force_min_allowed_freq(
        residue_freqs: Dict[NerpaResidue, Prob],
        min_allowed_freq: float
) -> Dict[NerpaResidue, Prob]:
    items: List[Tuple[NerpaResidue, Prob]] = (
        sorted(
            ((residue, prob)
            for residue, prob in residue_freqs.items()),
            key=lambda pair: pair[1]
        )
    )

    new_freqs: Dict[NerpaResidue, Prob] = {}
    missing_mass: float = 0.0
    for i, (residue, prob) in enumerate(items):
        tax = missing_mass / (len(items) - i)
        if prob - tax >= min_allowed_freq:
            for residue, prob in items[i:]:
                new_freqs[residue] = prob - tax
            break

        missing_mass += (min_allowed_freq - prob)
        new_freqs[residue] = min_allowed_freq
    else:
        raise ValueError(
            f"Could not force minimum allowed frequency of {min_allowed_freq} "
            f"on residue frequencies: {residue_freqs}"
        )

    return new_freqs
            
        

class NorineMonomersInfo:
    df: pl.DataFrame

    @classmethod
    def verify_schema(cls, df: pl.DataFrame) -> bool:
        return all(
            field in df.columns
            for field in NorineMonomersInfoSchema
        )

    def __init__(self, df: pl.DataFrame):
        if not NorineMonomersInfo.verify_schema(df):
            raise ValueError("DataFrame does not match the expected schema.")
        self.df = df

    @classmethod
    def from_tsv(cls, tsv_file: Path) -> NorineMonomersInfo:
        df = pl.read_csv(tsv_file, separator="\t")
        if not NorineMonomersInfo.verify_schema(df):
            raise ValueError("DataFrame does not match the expected schema.")
        return cls(df)

    def modification_frequencies(self) -> Dict[MODIFICATION_TYPE, Prob]:
        NS = NorineMonomersInfoSchema
        total_count: int = int(self.df.get_column(NS.COUNT).sum())
        methylated_count: int = int(
            self.df
            .filter(pl.col(NS.METHYLATED) == True)
            .get_column(NS.COUNT)
            .sum()
        )
        meth_prob: float = methylated_count / total_count

        total_known_chr_count: int = int(
            self.df
            .filter(pl.col(NS.CHIRALITY) != Chirality.UNKNOWN.name)
            .get_column(NS.COUNT)
            .sum()
        )
        d_chr_count: int = int(
            self.df
            .filter(pl.col(NS.CHIRALITY) == Chirality.D.name)
            .get_column(NS.COUNT)
            .sum()
        )
        d_chr_prob: float = d_chr_count / total_known_chr_count
            
        return {
            MODIFICATION_TYPE.METHYLATION: meth_prob,
            MODIFICATION_TYPE.EPIMERIZATION: d_chr_prob,
        }

    def residue_frequencies(
            self,
            supported_residues: set[NerpaResidue] | None = None,
            min_allowed_freq: float = 0.0
    ) -> Dict[NerpaResidue, Prob]:
        NS = NorineMonomersInfoSchema

        if supported_residues is None:
            supported_residues = set(
                self.df.get_column(NS.NERPA_RESIDUE)
            )

        total_count: int = int(self.df.get_column(NS.COUNT).sum())
        residue_freqs_df = (
            self.df
            .group_by(NS.NERPA_RESIDUE)
            .agg(pl.sum(NS.COUNT).alias("residue_count"))
            .with_columns(
                (pl.col("residue_count") / total_count)
                .alias("residue_prob")
            )
        )

        result: Dict[NerpaResidue, Prob] = {
            UNKNOWN_RESIDUE: 0.0
        }
        for row in residue_freqs_df.iter_rows(named=True):
            residue = NerpaResidue(row[NS.NERPA_RESIDUE])
            residue_prob = float(row["residue_prob"])
            if residue in supported_residues:
                result[residue] = residue_prob
            else:
                result[UNKNOWN_RESIDUE] += residue_prob

        for res in supported_residues:
            if res not in result:
                result[res] = 0.0

        result = force_min_allowed_freq(result, min_allowed_freq)
        return result           

