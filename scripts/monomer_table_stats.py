from pathlib import Path
import polars as pl
import yaml

nerpa_dir = Path(__file__).parent.parent
names_table_path = nerpa_dir / "data/nerpa_internal/monomer_names_table.tsv"
names_table = pl.read_csv(names_table_path, separator='\t')
# q: deduplicate names_table by MonomerName, and NameFormat, keeping the first occurrence
names_table = names_table.unique(subset=['MonomerName', 'NameFormat'])

# header: MonomerName     NameFormat      NerpaResidue    Modifications   Type    Comment

monomer_names_cfg_path = nerpa_dir / 'configs/monomers_config.yaml'

with open(monomer_names_cfg_path, 'r') as f:
    monomer_names_cfg = yaml.safe_load(f)

supported_residues = monomer_names_cfg['supported_residues']

# q: count the number of rows with NameFormat = 'rBAN/Norine' and NerpaResidue in supported_residues
rBAN_Norine_supported = names_table.filter(
    (pl.col('NameFormat') == 'rBAN/Norine') & 
    (pl.col('NerpaResidue').is_in(supported_residues))
)
print(rBAN_Norine_supported.head())
print(f"Number of rows with NameFormat = 'rBAN/Norine' and NerpaResidue in supported_residues: {rBAN_Norine_supported.height}")
