import sys
from pathlib import Path
import subprocess
import time
import webbrowser
import polars as pl

import yaml

def nrp_idx_to_id(nrp_idx: str) -> str:
    if any(ch.isalpha() for ch in nrp_idx):
              return nrp_idx
    total_len = 11 + len(nrp_idx.split('.')[1])
    zeros = '0' * (total_len - len('BGC') - len(nrp_idx))
    return 'BGC' + zeros + nrp_idx

def bgc_idx_to_id(bgc_idx: str) -> str:
    if any(ch.isalpha() for ch in bgc_idx):
              return bgc_idx
    total_len = 10
    zeros = '0' * (total_len - len('BGC') - len(bgc_idx))
    return 'BGC' + zeros + bgc_idx

def index_to_id(index: str) -> str:
    if any(ch.isalpha() for ch in index):
              return index
    if '.' in index:
        return nrp_idx_to_id(index)
    else:
        return bgc_idx_to_id(index)

args = sys.argv[1:]
open_html = '--dont-open-html' not in args
args = [arg for arg in args if arg != '--dont-open-html']

if len(args) == 1:
    nrp_id = index_to_id(args[0])
    bgc_id = nrp_id.split('.')[0]
elif len(args) == 2:
    bgc_id = index_to_id(args[0])
    nrp_id = index_to_id(args[1])
else:
    raise ValueError("Expected 1 or 2 arguments: <nrp_id> or <bgc_id> <nrp_id>")

print(f'Matching BGC {bgc_id} against NRP {nrp_id}')

nerpa_dir = Path(__file__).resolve().parent.parent
assert nerpa_dir.name.startswith('nerpa'), f"Expected nerpa_dir to be the nerpa repo, got {nerpa_dir}"
nerpa_script = nerpa_dir / "nerpa.py"
output_dir = nerpa_dir / "nerpa_results/trash_single"
output_dir.mkdir(parents=True, exist_ok=True)

local_paths = yaml.safe_load((nerpa_dir / "local_paths.yaml").read_text())
antismash_results_all = Path(local_paths['as_results_mibig4_nrps'])
antismash_results = antismash_results_all / bgc_id

pnrpdb = pl.read_csv(nerpa_dir / 'data/input/pnrpdb2_expanded.tsv', separator='\t')
smiles_tsv_path: Path = nerpa_dir / 'tmp' / f'{nrp_id}.tsv'
smiles_tsv_path.parent.mkdir(parents=True, exist_ok=True)
row_tbl: pl.DataFrame = pnrpdb.filter(pl.col('ID') == nrp_id).limit(1)
row_tbl.write_csv(smiles_tsv_path, separator='\t')

# Construct the command
command = [
    "python3", str(nerpa_script),
    "--antismash", str(antismash_results),
    "--smiles-tsv", smiles_tsv_path,
    "--col-id", "ID",
    "--output-dir", str(output_dir),
    "--force-output-dir",
    "--fast-matching",
    "--let-it-crash",
    "--dump-all-preprocessed",
    "--keep-intermediate-files",
    #"--draw-hmms"
]

subprocess.run(command, check=True)

if open_html:
    report_path = output_dir / 'report.html'
    # wait up to 30 seconds for the report to be created
    timeout = 30.0
    interval = 0.5
    elapsed = 0.0
    while elapsed < timeout and not report_path.exists():
        time.sleep(interval)
        elapsed += interval

    if report_path.exists():
        webbrowser.open(report_path.as_uri())
    else:
        print(f"report not found at {report_path}; skipping browser open")
