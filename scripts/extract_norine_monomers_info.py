from pathlib import Path
from typing import Literal, Tuple

import pandas as pd
from collections import Counter
import yaml

from src.monomer_names_helper import MonomerNamesHelper, NRP_Monomer, NRP_Monomer_Modification, UNKNOWN_MONOMER, \
    PKS_MONOMER, MonomerResidue, UNKNOWN_RESIDUE
from src.data_types import Chirality


def parsed_modifications(mods_str: str) -> Tuple[NRP_Monomer_Modification, ...]:
    def parse_mod(mod: str) -> NRP_Monomer_Modification:
        match mod:
            case 'MT':
                return NRP_Monomer_Modification.METHYLATION
            case other:
                return NRP_Monomer_Modification.UNKNOWN

    mods_list = eval(mods_str)
    return () if mods_list == [''] else tuple(map(parse_mod, mods_list))


def parsed_name(name: str,
                name_format: Literal['antismash', 'norine'],
                names_table: pd.DataFrame) -> NRP_Monomer:
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

    if name_format == 'norine' and name.startswith('X'):
        return UNKNOWN_MONOMER

    column_name = 'as_short' if name_format == 'antismash' else 'norine_rban_code'

    try:
        row = names_table[names_table[column_name] == name].iloc[0]
        residue = MonomerResidue(row['core'])
    except IndexError:
        return UNKNOWN_MONOMER  # TODO: raise ValueError when the table is ready
        # raise ValueError(f'Name {name} not found in the names table')

    return NRP_Monomer(residue=residue,
                       methylated=NRP_Monomer_Modification.METHYLATION in parsed_modifications(row['modifications']),
                       chirality=Chirality.UNKNOWN)



nerpa_dir = Path(__file__).parent.parent
# Load the Norine monomers table
norine_graphs_path = nerpa_dir / 'data/norine_graphs.csv'
norine_graphs_table = pd.read_csv(norine_graphs_path,
                                  delimiter=';')
graph_records = norine_graphs_table['graph'].to_list()

monomers_str = []
for graph_record in graph_records:
    monomers_str.extend(graph_record.split('@')[0].split(','))

monomers_with_chiralities = [(monomer[2:], Chirality.D) if monomer.startswith('D-') else (monomer, Chirality.L)
                             for monomer in monomers_str]

# Load monomer names helper
monomer_names_table_path = nerpa_dir / 'data/monomer_names_table.tsv'
monomer_names_table = pd.read_csv(monomer_names_table_path, sep='\t')

# Extract monomer names and modifications
parsed_monomers_with_ch = [(parsed_name(monomer, 'norine', monomer_names_table), ch)
                           for monomer, ch in monomers_with_chiralities]

# Count monomer residues
residues_counter = Counter(mon.residue for mon, ch in parsed_monomers_with_ch)
residues_total = sum(residues_counter.values())
residue_frequencies = {res: cnt / residues_total
                        for res, cnt in residues_counter.most_common()}

# output statistics
out_file = nerpa_dir / 'data/norine_monomers_info.yaml'
with out_file.open('w') as f:
    yaml.dump({
        'total_monomers': len(parsed_monomers_with_ch),
        'residue_frequencies': residue_frequencies,
        'methylated': sum(mon.methylated for mon, ch in parsed_monomers_with_ch),
        'epimerized': sum(ch == Chirality.D for mon, ch in parsed_monomers_with_ch)
    }, f)

