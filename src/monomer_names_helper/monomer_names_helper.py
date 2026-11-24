from __future__ import annotations
import ast
import math
from collections import defaultdict
from pathlib import Path
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
from src.general_type_aliases import Prob
from src.monomer_names_helper.monomer_names_table import (
    MonomerNameFormat,
    MonomerType,
    MonomerTag,
    MonomerNamesTableSchema,
    MonomerNamesTable_df,
    NerpaResidue,
    antiSMASH_MonomerName,
    NorineMonomerName,
    PARAS_MonomerName
)
import pandas as pd
import polars as pl
from dataclasses import dataclass
from enum import Enum, StrEnum, auto

import yaml

from src.pipeline.logging.logger import NerpaLogger


# Special "residues"
UNKNOWN_RESIDUE = NerpaResidue('unk')
NOT_NRPS_RESIDUE = NerpaResidue('not-nrps')  # residue that can't be attracted by an A domain, e.g. PKS
PKS_RESIDUE = NerpaResidue('pks')  # PKS residue in a PKS-NRPS hybrid
SPECIAL_RESIDUES = [UNKNOWN_RESIDUE, NOT_NRPS_RESIDUE, PKS_RESIDUE]

# Encoding of a monomer (including methylation, chirality, etc)
MonCode = NewType('MonCode', int)

class Chirality(Enum):
    L = auto()
    D = auto()
    UNKNOWN = auto()

    def __lt__(self, other):
        if not isinstance(other, Chirality):
            return NotImplemented
        order = {Chirality.L: 0, Chirality.D: 1, Chirality.UNKNOWN: 2}
        return order[self] < order[other]


class NRP_Monomer_Modification(Enum):  # post-translational modification
    METHYLATION = auto()
    UNKNOWN = auto()

yaml.add_representer(Chirality, enum_representer)
yaml.add_representer(NRP_Monomer_Modification, enum_representer)

@dataclass
class NRP_Monomer:
    residue: NerpaResidue
    methylated: bool = False
    chirality: Chirality = Chirality.UNKNOWN
    is_pks_hybrid: bool = False

    def __hash__(self):
        return hash((self.residue, self.methylated, self.chirality, self.is_pks_hybrid))

    def __lt__(self, other):
        if other is None:
            return False
        return (self.residue, self.methylated, self.chirality.name, self.is_pks_hybrid) < \
               (other.residue, other.methylated, other.chirality.name, other.is_pks_hybrid)

    def to_str_short(self) -> str:
        return (
            f"{self.residue}"
            + ("_Me" if self.methylated else "")
            + f"_Chr-{self.chirality.name}"
            + ("_PKS" if self.is_pks_hybrid else "")
        )

    @classmethod
    def from_str_short(cls, s: str) -> NRP_Monomer:
        parts = s.split('_')
        # quick fix to handle special residues that start with '_'
        if s.startswith('_'):
            residue = NerpaResidue('_' + parts[1])
            rest = parts[2:]
        else:
            residue = NerpaResidue(parts[0])
            rest = parts[1:]
            
        methylated = any(part == "Me" for part in rest)
        chirality = next(
            (Chirality[part.split('-')[1]]
             for part in rest
             if part.startswith("Chr-")),
            Chirality.UNKNOWN
        )
        is_pks_hybrid = any(part == "PKS" for part in rest)
        return cls(
            residue=residue,
            methylated=methylated,
            chirality=chirality,
            is_pks_hybrid=is_pks_hybrid
        )


UNKNOWN_MONOMER = NRP_Monomer(
    residue=UNKNOWN_RESIDUE,
    methylated=False,  # TODO: should we have a separate UNKNOWN_METHYLATION?,
    chirality=Chirality.UNKNOWN
)
PKS_MONOMER = NRP_Monomer(
    residue=PKS_RESIDUE,
    methylated=False,  # TODO: should we have a separate UNKNOWN_METHYLATION?,
    chirality=Chirality.UNKNOWN
)
NOT_NRPS_MONOMER = NRP_Monomer(
    residue=NOT_NRPS_RESIDUE,
    methylated=False,  # TODO: should we have a separate UNKNOWN_METHYLATION?,
    chirality=Chirality.UNKNOWN
)


class MonomersDefaultFrequencies(NamedTuple):
    residue: Dict[NerpaResidue, Prob]
    methylation: Prob
    d_chirality: Prob


class MonomerNamesHelper:
    names_table: MonomerNamesTable_df
    supported_residues: List[NerpaResidue]

    _parsed_name_cache: Dict[Tuple[str, str], NRP_Monomer]

    mon_to_int: Dict[NRP_Monomer, MonCode]
    int_to_mon: Dict[MonCode, NRP_Monomer]


    def __init__(
            self,
            names_table_path: Path,
            supported_residues: Optional[List[NerpaResidue]] = None
    ):
        SC = MonomerNamesTableSchema  # for brevity

        self._parsed_name_cache = {}
        self.names_table = MonomerNamesTable_df.from_tsv(names_table_path)
        self.supported_residues = (
            supported_residues
            if supported_residues is not None
            else [
                NerpaResidue(residue)
                for residue in self.names_table.df[SC.NERPA_RESIDUE].unique()
            ]
        )

        # Check that all supported residues are present in the names table.
        table_residues = set(
            self.names_table.df[SC.NERPA_RESIDUE]
        )
        assert all(res in table_residues for res in self.supported_residues), \
            f"Some supported residues are absent in the names table: " \
            f"{set(self.supported_residues) - table_residues}"

        # Compute monomer integer codes and mappings
        self.int_to_mon = {}
        self.mon_to_int = {}
        for mon_res_int, mon_res in enumerate(self.supported_residues + SPECIAL_RESIDUES):
            for meth_int, methylated in enumerate([False, True]):
                for chir_int, chirality in enumerate([Chirality.L, Chirality.D, Chirality.UNKNOWN]):
                    for is_pks_hybrid_int, is_pks_hybrid in enumerate([False, True]):
                        mon = NRP_Monomer(
                            residue=NerpaResidue(mon_res),
                            methylated=methylated,
                            chirality=chirality,
                            is_pks_hybrid=is_pks_hybrid
                        )
                        mon_int = mon_res_int * 12 + meth_int * 6 + chir_int * 2 + is_pks_hybrid_int
                        self.mon_to_int[mon] = MonCode(mon_int)
                        self.int_to_mon[MonCode(mon_int)] = mon

    def is_proper_monomer(self, mon: NRP_Monomer) -> bool:
        return all([
            not mon.is_pks_hybrid,
            mon.chirality != Chirality.UNKNOWN,
            mon.residue not in [PKS_RESIDUE, NOT_NRPS_RESIDUE]
        ])

    def proper_monomers(self) -> Iterable[NRP_Monomer]:
        """
        Get a list of all proper monomers (excluding UNKNOWN and PKS).
        """
        return filter(self.is_proper_monomer, self.mon_to_int.keys())

    def get_table_entry(
            self,
            name: str, 
            name_format: MonomerNameFormat
    ) -> Optional[dict[MonomerNamesTableSchema, Any]]:
        SC = MonomerNamesTableSchema  # for brevity
        rows = (
            self.names_table.df
            .filter((pl.col(SC.MONOMER_NAME) == name)
                    & (pl.col(SC.NAME_FORMAT) == name_format))
        )

        if rows.is_empty():
            return None

        assert rows.height == 1, f"Expected exactly one row for {name} in {name_format}, but got {rows.height}"
        row = rows.to_dicts()[0]
        return {
            SC(col): SC(col).to_python_type(row[col])
            for col in row
        }


    def parsed_name(
            self, name: str,
            name_format: MonomerNameFormat,
            log: Optional[NerpaLogger] = None
    ) -> NRP_Monomer:
        SC = MonomerNamesTableSchema  # for brevity
        if (name, name_format) in self._parsed_name_cache:
            return self._parsed_name_cache[(name, name_format)]

        if name_format == MonomerNameFormat.RBAN_NORINE and name.startswith('X'):
            return UNKNOWN_MONOMER

        row = self.get_table_entry(name, name_format)
        if row is None:
            return UNKNOWN_MONOMER

        mon_type = row[SC.TYPE]
        if mon_type == MonomerType.PKS:
            return PKS_MONOMER

        # if mon_type == 'NOT_NRPS':  # Dil and Dap were erroneously marked as NOT_NRPS. As a quick fix -- treat all monomers as NRPS
        #     return NOT_NRPS_MONOMER


        residue = row[SC.NERPA_RESIDUE]
        if residue not in self.supported_residues:
            residue = UNKNOWN_RESIDUE
        try:
            mods = row[SC.MODIFICATIONS]
            methylated = any(
                mod_str in mods
                for mod_str in ['MT', 'OMe', 'NMe', 'MeO']
            )
            chirality = Chirality.D if 'E' in mods else Chirality.L   # TODO: make a separate column for chirality
        except:
            print(f'Error parsing modifications for {name}\n'
                  f'Table row: {row}\n')
            raise  # debug

        self._parsed_name_cache[(name, name_format)] = NRP_Monomer(
            residue=residue,
            methylated=methylated,
            chirality=chirality,
            is_pks_hybrid=(mon_type == MonomerType.PKS_HYBRID)
        )
        return self._parsed_name_cache[(name, name_format)]

    def all_residues(self) -> List[NerpaResidue]:
        return self.supported_residues + SPECIAL_RESIDUES

