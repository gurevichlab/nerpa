import sys
from pathlib import Path
import subprocess
import time
import webbrowser
import polars as pl
import argparse

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

def parse_args(nerpa_dir: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="run Nerpa on an MIBiG BGC vs its NRP(s)")
    parser.add_argument("bgc_id", type=str, help="BGC ID to process")
    parser.add_argument(
        "--dont-open-html",
        action="store_true",
        help="Do not open the generated report.html in a browser",
    )
    parser.add_argument(
        "-o", "--output-dir",
        dest="output_dir",
        type=Path,
        default=nerpa_dir / "nerpa_results/mibig_bgc_vs_its_nrps",
        help="Output directory (default: nerpa_results/mibig_bgc_vs_its_nrps)"
    )

    return parser.parse_args()

    
def main():
    nerpa_dir = Path(__file__).resolve().parent.parent
    assert nerpa_dir.name.startswith('nerpa'), f"Expected nerpa_dir to be the nerpa repo, got {nerpa_dir}"
    nerpa_script = nerpa_dir / "nerpa.py"

    args = parse_args(nerpa_dir)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    bgc_id = index_to_id(args.bgc_id)

    pnrpdb = pl.read_csv(nerpa_dir / 'data/input/pnrpdb2_expanded.tsv', separator='\t')
    smiles_tsv_path: Path = nerpa_dir / 'tmp' / f'{bgc_id}_nrps.tsv'
    smiles_tsv_path.parent.mkdir(parents=True, exist_ok=True)
    nrps_table: pl.DataFrame = pnrpdb.filter(pl.col('ID').str.starts_with(bgc_id))
    nrps_table.write_csv(smiles_tsv_path, separator='\t')

    nrps = nrps_table['ID'].to_list()

    print(f'Matching BGC {bgc_id} against {len(nrps)} NRPs:\n{nrps}')


    local_paths = yaml.safe_load((nerpa_dir / "local_paths.yaml").read_text())
    antismash_results_all = Path(local_paths['as_results_mibig4_nrps'])
    antismash_results = antismash_results_all / bgc_id

    # Construct the command
    command = [
        "python3", str(nerpa_script),
        "--antismash", str(antismash_results),
        "--smiles-tsv", smiles_tsv_path,
        "--col-id", "ID",
        "--output-dir", str(args.output_dir),
        "--force-output-dir",
        "--let-it-crash",
        "--dump-all-preprocessed",
        "--draw-hmms"
    ]

    subprocess.run(command, check=True)

    if not args.dont_open_html:
        report_path = args.output_dir / 'report.html'
        # wait up to 30 seconds for the report to be created
        timeout = 30.0
        interval = 0.5
        elapsed = 0.0
        while elapsed < timeout and not report_path.exists():
            time.sleep(interval)
            elapsed += interval

        if report_path.exists():
            webbrowser.open(report_path.resolve().as_uri())
        else:
            print(f"report not found at {report_path}; skipping browser open")

if __name__ == '__main__':
    main()
