import math
from collections import defaultdict
from math import log
from typing import (
    Any,
    Dict,
    Literal,
    List,
    Optional,
    NamedTuple,
    NewType,
    Tuple, Iterable
)
from src.build_output.yaml_representers import enum_representer
from src.general_type_aliases import Prob, LogProb
import pandas as pd
import polars as pl
from dataclasses import dataclass
from enum import Enum, auto
from src.monomer_names_helper.monomer_names_helper import (
    MonomerNamesHelper,
    NRP_Monomer,
    NerpaResidue,
    Chirality,
    UNKNOWN_RESIDUE,
    SPECIAL_RESIDUES,
)
import polars as pl

class MonomerFrequenciesHelper:
    residue_freqs: Dict[NerpaResidue, Prob]
    methylation_freq: Prob
    d_chirality_freq: Prob

    supported_residues: List[NerpaResidue]

    def __init__(
            self,
            norine_monomers_info: pl.DataFrame,
            supported_residues: Optional[List[NerpaResidue]],
            missing_residue_freq: Optional[Prob] = None
    ):
        if missing_residue_freq is None:
            missing_residue_freq = Prob(0.0)  # Default to 0 if not provided
        if supported_residues is None:
            supported_residues = [
                row['Nerpa_residue']
                for row in norine_monomers_info.iter_rows(named=True)
            ]
        
        self.supported_residues = supported_residues

        # 1. Compute residue_freqs
        residue_cnts: Dict[NerpaResidue, int] = defaultdict(int)
        for row in norine_monomers_info.iter_rows(named=True):
            res = (
                row['Nerpa_residue']
                if row['Nerpa_residue'] in supported_residues
                else UNKNOWN_RESIDUE
            )
            residue_cnts[res] += row['count']

        missing_residues = [
            res
            for res in supported_residues
            if res not in residue_cnts
        ]
        total_weight_of_missing = missing_residue_freq * len(missing_residues)

        total_monomers = sum(residue_cnts.values())
        self.residue_freqs = {
            res: (
                Prob((residue_cnts[res] / total_monomers) * (1 - total_weight_of_missing))
                if res in residue_cnts
                else missing_residue_freq
            )
            for res in supported_residues + SPECIAL_RESIDUES
        }

        # 2. Modification frequencies
        total_methylated = sum(
            row['count'] for row in norine_monomers_info.iter_rows(named=True)
            if row['methylated']
        )
        self.methylation_freq = Prob(total_methylated / total_monomers)

        total_d_chirality = sum(
            row['count'] for row in norine_monomers_info.iter_rows(named=True)
            if row['chirality'] == Chirality.D.name
        )
        self.d_chirality_freq = Prob(total_d_chirality / total_monomers)

        
    def monomer_default_freq(
            self,
            mon: NRP_Monomer,
            handle_pks_hybrids: bool = False
    ) -> Prob:
        return math.e ** self.monomer_default_log_freq(mon, handle_pks_hybrids)

    def monomer_default_log_freq(
            self,
            mon: NRP_Monomer,
            handle_pks_hybrids: bool = False
    ) -> LogProb:
        res = (
            mon.residue
            if not mon.is_pks_hybrid or handle_pks_hybrids
            else UNKNOWN_RESIDUE
        )
        residue_log_freq = log(self.residue_freqs[res])
        methylation_log_freq = (
            log(self.methylation_freq)
            if mon.methylated
            else log(1 - self.methylation_freq)
        )

        d_chr_freq = self.d_chirality_freq
        match mon.chirality:
            case Chirality.D:
                chirality_log_freq = log(d_chr_freq)
            case Chirality.L:
                chirality_log_freq = log(1 - d_chr_freq)
            case Chirality.UNKNOWN:
                chirality_log_freq = (
                    d_chr_freq * log(d_chr_freq)
                    + (1 - d_chr_freq) * log(1 - d_chr_freq)
                )

        return residue_log_freq + methylation_log_freq + chirality_log_freq
