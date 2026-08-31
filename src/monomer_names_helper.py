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
from src.general_type_aliases import Prob
import pandas as pd
import polars as pl
from dataclasses import dataclass
from enum import Enum, auto

import yaml

from src.pipeline.logging.logger import NerpaLogger

antiSMASH_MonomerName = NewType('antiSMASH_MonomerName', str)  # antiSMASH Short
NorineMonomerName = NewType('NorineMonomerName', str)
PARAS_MonomerName = NewType('PARAS_MonomerName', str)
NerpaResidue = NewType('NerpaResidue', str)

# Special "residues"
UNKNOWN_RESIDUE = NerpaResidue('_UNKNOWN')
NOT_NRPS_RESIDUE = NerpaResidue('_NOT_NRPS')  # residue that can't be attracted by an A domain, e.g. PKS
PKS_RESIDUE = NerpaResidue('_PKS')  # PKS residue in a PKS-NRPS hybrid
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
            f'{self.residue}'
            "_Me" if self.methylated else ""
            f"_Chr-{self.chirality.name}"
            "_PKS" if self.is_pks_hybrid else ""
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
    '''
    names_table: a dataframe with columns
     MonomerName,NameFormat,antiSMASH_short_core/Nerpa_residue,Modifications,Type,Comment
    '''
    names_table: pl.DataFrame
    supported_residues: List[NerpaResidue]

    _parsed_name_cache: Dict[Tuple[str, str], NRP_Monomer]

    mon_to_int: Dict[NRP_Monomer, MonCode]
    int_to_mon: Dict[MonCode, NRP_Monomer]


    def __init__(
            self,
            names_table: pl.DataFrame,
            supported_residues: Optional[List[NerpaResidue]] = None
    ):
        self._parsed_name_cache = {}
        self.names_table = names_table
        self.supported_residues = (
            supported_residues
            if supported_residues is not None
            else list(self.names_table['NerpaResidue'].unique())
        )

        # Check that all supported residues are present in the names table
        assert all(res in self.names_table['NerpaResidue']
                   for res in self.supported_residues), \
            f"Some supported residues are absent in the names table: " \
            f"{set(self.supported_residues) - set(self.names_table['NerpaResidue'])}"

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

    def parsed_name(self, name: str,
                    name_format: Literal['antiSMASH_short', 'rBAN/Norine', 'PARAS'],
                    log: Optional[NerpaLogger] = None) -> NRP_Monomer:
        if name_format not in ('antiSMASH_short', 'rBAN/Norine', 'PARAS'):
            raise ValueError(f'Unsupported name format: {name_format}')
        if (name, name_format) in self._parsed_name_cache:
            return self._parsed_name_cache[(name, name_format)]

        if name_format == 'rBAN/Norine' and name.startswith('X'):
            return UNKNOWN_MONOMER

        rows = (
            self.names_table
            .filter((pl.col('MonomerName') == name)
                    & (pl.col('NameFormat') == name_format))
        )

        if rows.is_empty():  # monomer name not found
            # print(f'WARNING: Entry {(name, name_format)} not found in the names table. Parsing as UNKNOWN_MONOMER')
            '''
            if log is not None:
                log.error(f'Name {name} not found in the names table. Parsing as UNKNOWN_MONOMER')
            else:
                #with open('missing_norine_names.txt', 'a') as f:
                #    f.write(f'{name}\n')
                print(f'WARNING: Name {name} not found in the names table. Parsing as UNKNOWN_MONOMER')
            '''
            return UNKNOWN_MONOMER

        row = rows.to_dicts()[0]  # it should be exactly one row

        mon_type = row['Type']
        if mon_type == 'PKS':
            return PKS_MONOMER

        # if mon_type == 'NOT_NRPS':  # Dil and Dap were erroneously marked as NOT_NRPS. As a quick fix -- treat all monomers as NRPS
        #     return NOT_NRPS_MONOMER

        assert mon_type in ('NRPS', 'NA', 'PKS_Hybrid', 'NOT_NRPS'), \
            f'Unexpected monomer type: {name} -- {mon_type}'

        residue = NerpaResidue(row['NerpaResidue'])
        if residue not in self.supported_residues:
            residue = UNKNOWN_RESIDUE
        try:
            mods = eval(row['Modifications'])
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
            is_pks_hybrid=mon_type == 'PKS_Hybrid'
        )
        return self._parsed_name_cache[(name, name_format)]

    def all_residues(self) -> List[NerpaResidue]:
        return self.supported_residues + SPECIAL_RESIDUES

