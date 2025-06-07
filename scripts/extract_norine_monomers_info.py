from pathlib import Path

import pandas as pd
from collections import Counter
import yaml

nerpa_dir = Path(__file__).parent.parent
# Load the Norine monomers table
norine_graphs_path = nerpa_dir / 'data/norine_graphs.csv'
norine_graphs_table = pd.read_csv(norine_graphs_path,
                                  delimiter=';')
graph_records = norine_graphs_table['graph'].to_list()

monomers_str = []
for graph_record in graph_records:
    monomers_str.extend(graph_record.split('@')[0].split(','))

# Count monomer cnts
monomers_counter = Counter(monomers_str)

# output statistics
out_file = nerpa_dir / 'data/norine_monomers_cnts.yaml'
with out_file.open('w') as f:
    yaml.dump(dict(monomers_counter), f)

