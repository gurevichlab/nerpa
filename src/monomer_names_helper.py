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


antiSMASH_MonomerName = NewType('antiSMASH_MonomerName', str)  # antiSMASH Short
NorineMonomerName = NewType('NorineMonomerName', str)
MonomerResidue = NewType('MonomerResidue', str)
UNKNOWN_RESIDUE = MonomerResidue('unknown')


class NRP_Monomer_Modification(Enum):  # post-translational modification
    METHYLATION = auto()
    UNKNOWN = auto()


class ParsedMonomerName(NamedTuple):
    residue: MonomerResidue
    modifications: Tuple[NRP_Monomer_Modification, ...]


@dataclass
class MonomerNamesHelper:
    names_table: pd.DataFrame

    def my_parser(self, name: str) -> ParsedMonomerName:  # stub before monomer table is finished
        *mods, res = name.split('-')
        if len(mods) == 0:
            modifications = ()
        else:
            if any(mod in ['Me', 'CMe', 'NMe', 'OMe'] for mod in mods):
                modifications = (NRP_Monomer_Modification.METHYLATION,)
            else:
                modifications = (NRP_Monomer_Modification.UNKNOWN,)

        if res == 'aThr/Thr':
            res = 'Thr'
        if res == 'Ile/aIle':
            res = 'Ile'
        if res not in self.names_table['core'].values:
            res = UNKNOWN_RESIDUE
        return ParsedMonomerName(residue=MonomerResidue(res),
                                 modifications=modifications)


    @classmethod
    def parsed_modifications(cls, mods_str: str) -> Tuple[NRP_Monomer_Modification, ...]:
        def parse_mod(mod: str) -> NRP_Monomer_Modification:
            match mod:
                case 'MT': return NRP_Monomer_Modification.METHYLATION
                case other: return NRP_Monomer_Modification.UNKNOWN

        mods_list = eval(mods_str)
        return () if mods_list == ['-'] else tuple(map(parse_mod, mods_list))

    def parsed_name(self, name: str, name_format: Literal['antismash', 'norine']) -> ParsedMonomerName:
        if name_format == 'norine' and name.startswith('X'):
            return ParsedMonomerName(residue=UNKNOWN_RESIDUE, modifications=())

        column_name = 'as_short' if name_format == 'antismash' else 'norine_rban_code'

        try:
            row = self.names_table[self.names_table[column_name] == name].iloc[0]
        except IndexError:
            return ParsedMonomerName(residue=UNKNOWN_RESIDUE, modifications=())  # TODO: raise ValueError when the table is ready
            # raise ValueError(f'Name {name} not found in the names table')

        return ParsedMonomerName(residue=MonomerResidue(row['core']),
                                 modifications=self.parsed_modifications(row['modifications']))
