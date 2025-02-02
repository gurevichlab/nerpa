from typing import (
    Any,
    Dict,
    Literal,
    List,
    Optional,
    NamedTuple,
    NewType,
    Tuple
)
import pandas as pd
from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path
import yaml


antiSMASH_MonomerName = NewType('antiSMASH_MonomerName', str)  # antiSMASH Short
NorineMonomerName = NewType('NorineMonomerName', str)
MonomerResidue = NewType('MonomerResidue', str)
UNKNOWN_RESIDUE = MonomerResidue('unknown')
PKS = MonomerResidue('PKS')

def enum_representer(dumper, e: Enum):
    return dumper.represent_scalar(f'!{e.__class__.__name__}', e.name)


class Chirality(Enum):
    L = auto()
    D = auto()
    UNKNOWN = auto()


class NRP_Monomer_Modification(Enum):  # post-translational modification
    METHYLATION = auto()
    UNKNOWN = auto()

yaml.add_representer(Chirality, enum_representer)
yaml.add_representer(NRP_Monomer_Modification, enum_representer)

@dataclass
class NRP_Monomer:
    residue: MonomerResidue
    methylated: bool = False
    chirality: Chirality = Chirality.UNKNOWN
    is_pks_hybrid: bool = False

    def __hash__(self):
        return hash((self.residue, self.methylated, self.chirality, self.is_pks_hybrid))


@dataclass
class MonomerNamesHelper:
    names_table: pd.DataFrame
    supported_residues: List[MonomerResidue]
    mon_to_int: Dict[NRP_Monomer, int] = None
    int_to_mon: Dict[int, NRP_Monomer] = None
    pks_names: List[NorineMonomerName] = None  # placeholder for PKS names

    def __post_init__(self):
        if self.pks_names is None:
            self.pks_names = []

        self.int_to_mon = {}
        self.mon_to_int = {}
        for mon_res_int, mon_res in enumerate(self.supported_residues):
            for meth_int, methylated in enumerate([False, True]):
                for chir_int, chirality in enumerate([Chirality.L, Chirality.D, Chirality.UNKNOWN]):
                    for is_pks_hybrid_int, is_pks_hybrid in enumerate([False, True]):
                        mon = NRP_Monomer(residue=MonomerResidue(mon_res),
                                          methylated=methylated,
                                          chirality=chirality,
                                          is_pks_hybrid=is_pks_hybrid)
                        mon_int = mon_res_int * 12 + meth_int * 6 + chir_int * 2 + is_pks_hybrid_int
                        self.mon_to_int[mon] = mon_int
                        self.int_to_mon[mon_int] = mon


    @classmethod
    def parsed_modifications(cls, mods_str: str) -> Tuple[NRP_Monomer_Modification, ...]:
        def parse_mod(mod: str) -> NRP_Monomer_Modification:
            match mod:
                case 'MT': return NRP_Monomer_Modification.METHYLATION
                case other: return NRP_Monomer_Modification.UNKNOWN

        mods_list = eval(mods_str)
        return () if mods_list == [''] else tuple(map(parse_mod, mods_list))

    def parsed_name(self, name: str, name_format: Literal['antismash', 'norine']) -> NRP_Monomer:
        if name == 'Ile/aIle':
            name = 'Ile'  # temporary fix for the table
        if name == 'NMe-aIle/NMe-Ile':
            name = 'NMe-Ile'
        if name == 'aThr/Thr':
            name = 'Thr'
        if name == 'NMe-aThr/NMe-Thr':
            name = 'NMe-Thr'
        if name == 'NMe-bMe-Leu/diMe-aIle':
            name = 'NMe-Leu'

        UNKNOWN_MONOMER = NRP_Monomer(residue=UNKNOWN_RESIDUE,
                                      methylated=False,  # TODO: should we have a separate UNKNOWN_METYLATION?,
                                      chirality=Chirality.UNKNOWN)
        PKS_MONOMER = NRP_Monomer(residue=PKS,
                                  methylated=False,  # TODO: should we have a separate UNKNOWN_METYLATION?,
                                  chirality=Chirality.UNKNOWN)
        if name_format == 'norine' and name.startswith('X'):
            return UNKNOWN_MONOMER
        if name_format == 'norine' and name in self.pks_names:
            return PKS_MONOMER

        column_name = 'as_short' if name_format == 'antismash' else 'norine_rban_code'

        try:
            row = self.names_table[self.names_table[column_name] == name].iloc[0]
            residue = MonomerResidue(row['core']) if row['core'] in self.supported_residues else UNKNOWN_RESIDUE
        except IndexError:
            return UNKNOWN_MONOMER  # TODO: raise ValueError when the table is ready
            # raise ValueError(f'Name {name} not found in the names table')

        return NRP_Monomer(residue=residue,
                           methylated=NRP_Monomer_Modification.METHYLATION in self.parsed_modifications(row['modifications']),
                           chirality=Chirality.UNKNOWN)
