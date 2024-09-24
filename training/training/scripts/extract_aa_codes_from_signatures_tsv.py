from pathlib import Path
from collections import defaultdict
import yaml

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
        try:
            row = self.names_table[self.names_table[column_name] == name].iloc[0]
        except Exception:
            raise
        return ParsedMonomerName(residue=MonomerResidue(row['core']),
                                 modifications=self.parsed_modifications(row['modifications']))


monomer_names_tsv = Path('/home/ilianolhin/git/nerpa2/configs/monomers_tables_v2/monomers_short_unique.tsv')
monomer_names_df = pd.read_csv(monomer_names_tsv, sep='\t')
names_helper = MonomerNamesHelper(monomer_names_df)

sig_tsv = Path('/home/ilianolhin/git/nerpa2/configs/signatures.tsv')
sig_df = pd.read_csv(sig_tsv, sep='\t')

sig_df.columns = ['aa10', 'aa34', 'validated_cores_loose', 'validated_cores_strict', 'genomes_ids']
aa_codes = {'aa10': {}, 'aa34': {}}
for _, row in sig_df.iterrows():
    for core in row['validated_cores_loose'].split('|'):
        residue = names_helper.parsed_name(core, 'antismash').residue
        if residue in aa_codes['aa10']:
            if row['aa10'] not in aa_codes['aa10'][residue]:
                aa_codes['aa10'][residue].append(row['aa10'])
        else:
            aa_codes['aa10'][residue] = [row['aa10']]

        if residue in aa_codes['aa34']:
            if row['aa34'] not in aa_codes['aa34'][residue]:
                aa_codes['aa34'][residue].append(row['aa34'])
        else:
            aa_codes['aa34'][residue] = [row['aa34']]

def write_yaml(data, out_file: Path):
    # dirty hack to erase information about types and make output less verbose
    # https://github.com/yaml/pyyaml/issues/408
    yaml.emitter.Emitter.prepare_tag = lambda self, tag: ''

    # another hack (albeit less dirty) to forbid yaml creating references
    # https://stackoverflow.com/questions/13518819/avoid-references-in-pyyaml
    yaml.Dumper.ignore_aliases = lambda *args: True

    with open(out_file, 'w') as out:
        yaml.dump(data, out,
                  default_flow_style=None, sort_keys=False)

write_yaml(aa_codes, Path("../../../configs/aa_codes.yaml"))