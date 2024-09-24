from src.monomer_names_helper import MonomerNamesHelper
from pathlib import Path
import csv
from collections import Counter
import pandas as pd
from pathlib import Path

# this script is in nerpa_dir/src/scripts/calculate_norine_monomers_freqs.py
nerpa_dir = Path(__file__).parent.parent.parent
monomer_names_table = pd.read_csv(nerpa_dir / 'configs' / 'monomer_names_table.tsv',
                                  sep='\t')
monomer_names_helper = MonomerNamesHelper(monomer_names_table)

norine_table = pd.read_csv(nerpa_dir / 'data/norine_graphs.csv', sep=';')
# q: set column names to ['norine_id', 'norine_familty', 'graph_representation']
norine_table.columns = ['norine_id',  'graph_representation']
norine_table = norine_table.replace({float('nan'): None})

norine_monomers = []
for _, row in norine_table.iterrows():
    if row['graph_representation'] is None:
        continue
    norine_graph_str = row['graph_representation']
    monomers_str = norine_graph_str.split('@')[0]
    norine_monomers.extend(monomers_str.split(','))

norine_residues = [monomer_names_helper.parsed_name(monomer_name, 'norine').residue
                   for monomer_name in norine_monomers]
residues_counter = Counter(norine_residues)
# q: write residues frequencies fractions
residues_freqs = {residue: count / len(norine_residues) for residue, count in residues_counter.most_common()}
with open(nerpa_dir / 'data/norine_residues_freqs.tsv', 'w') as file:
    file.write('residue\tfrequency\n')
    for residue, freq in residues_freqs.items():
        file.write(f'{residue}\t{freq}\n')
