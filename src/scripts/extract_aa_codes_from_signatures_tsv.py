from src.write_results import write_yaml
from pathlib import Path
import pandas as pd
import yaml

sig_tsv = Path('/home/ilianolhin/git/nerpa2/configs/signatures.tsv')
sig_df = pd.read_csv(sig_tsv, sep='\t')

# q: annotate columns with the names ['aa10', 'aa34', 'validated_cores_loose', 'validated_cores_strict', 'genomes_ids']
# q: what is the type of the 'genomes_ids' column? All types are str

sig_df.columns = ['aa10', 'aa34', 'validated_cores_loose', 'validated_cores_strict', 'genomes_ids']
aa_codes = {'aa10': {}, 'aa34': {}}
for _, row in sig_df.iterrows():
    for core in row['validated_cores_loose'].split('|'):
        if core in aa_codes['aa10']:
            if row['aa10'] not in aa_codes['aa10'][core]:
                aa_codes['aa10'][core].append(row['aa10'])
        else:
            aa_codes['aa10'][core] = [row['aa10']]

        if core in aa_codes['aa34']:
            if row['aa34'] not in aa_codes['aa34'][core]:
                aa_codes['aa34'][core].append(row['aa34'])
        else:
            aa_codes['aa34'][core] = [row['aa34']]

write_yaml(aa_codes, Path('/home/ilianolhin/git/nerpa2/configs/aa_codes.yaml'))
