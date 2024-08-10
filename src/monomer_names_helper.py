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

        column_name = 'as_short' if name_format == 'antismash' else 'as_norine'
        row = self.names_table[self.names_table[column_name] == name]
        return ParsedMonomerName(residue=MonomerResidue(row['core']),
                                 modifications=self.parsed_modifications(row['modifications']))




