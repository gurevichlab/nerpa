import os
import sys

import pandas as pd
import json
import subprocess
from pathlib import Path
from typing import Dict
import argparse

import yaml

from src.rban_parsing.rban_parser import Parsed_rBAN_Record


def run_nerpa_results_on_pnrpdb2(nerpa_dir: Path,
                                 pnrpdb2_path: Path,
                                 output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    nerpa_script = nerpa_dir / "nerpa.py"

    num_threads = 32 if os.cpu_count() >= 40 else 8  # to work on cluster and locally
    # Construct the command
    command = [
        "python3", str(nerpa_script),
        "--antismash", str(nerpa_dir / "test_data" / "antismash"),
        "--smiles-tsv", str(pnrpdb2_path),
        "--col-id", "ID",
        "--output-dir", str(output_dir),
        "--force-output-dir",
        "--max-num-matches-per-nrp", '1',
        "--min-num-matches-per-nrp", '1',
        "--max-num-matches", '0',
        "--skip-molecule-drawing",
        "--threads", str(num_threads),
        "--fast-matching",
        "--dump-all-preprocessed",
        "--keep-intermediate-files",
    ]

    # Execute the command and capture output
    print('Executing command:', ' '.join(command))
    try:
        result = subprocess.run(command, check=True)
        print(f"Nerpa run completed successfully")
    except subprocess.CalledProcessError as e:
        print(f"Error running Nerpa:\n{e.stderr}")
        raise

def load_nrp_to_representative(nerpa_results: Path) -> Dict[str, str]:
    import yaml
    nrp_representatives_path = (nerpa_results
                                / 'preprocessed_input'
                                / 'NRP_representatives.yaml')
    with nrp_representatives_path.open() as f:
        nrp_representatives = yaml.safe_load(f)
    nrp_id_to_repr_id = {entry['nrp_id'][0]: entry['representative_id'][0]
                         for entry in nrp_representatives}
    return nrp_id_to_repr_id


def parse_args(nerpa_dir: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Preprocess pnrpdb2 for nerpa results and input tables.")
    parser.add_argument(
        "--use-existing-results",
        action='store_true',
        help="Skip running nerpa on pnrpdb2 and use existing results (for debugging)."
    )
    parser.add_argument(
        "--pnrpdb2-path",
        type=Path,
        help="Path to pnrpdb2 TSV file (if not provided, will use data/input/pnrpdb2.tsv)."
    )

    args = parser.parse_args()
    if args.pnrpdb2_path is None:
        args.pnrpdb2_path = (nerpa_dir
                             / 'data'
                             / 'input'
                             / 'pnrpdb2.tsv')
    return args

def main():
    nerpa_dir = Path(__file__).resolve().parent.parent
    args = parse_args(nerpa_dir)
    pnrpdb2_path = args.pnrpdb2_path
    pnrpdb_pref = pnrpdb2_path.stem
    print(f'Using pnrpdb2 path: {pnrpdb2_path}, prefix for output tables: {pnrpdb_pref}')
    pnrpdb2_df = pd.read_csv(pnrpdb2_path, sep='\t')

    nerpa_results = (nerpa_dir
                     / 'nerpa_results'
                     / 'preprocess_pnrpdb2')

    if not args.use_existing_results:
        run_nerpa_results_on_pnrpdb2(nerpa_dir=nerpa_dir,
                                     pnrpdb2_path=pnrpdb2_path,
                                     output_dir=nerpa_results)

    nrp_id_to_repr_id = load_nrp_to_representative(nerpa_results)
    rban_records = [
        Parsed_rBAN_Record.from_dict(record)
        for record in yaml.safe_load((nerpa_results / 'preprocessed_input' / 'parsed_rban_records.yaml').open())
    ]

    def is_mibig_norine(nrp_id: str) -> bool:
        return nrp_id.startswith('BGC') or nrp_id.startswith('NOR')

    def is_representative(nrp_id: str) -> bool:
        return nrp_id in nrp_id_to_repr_id.values()

    # 1. pnrpdb2_mibig_norine.tsv
    pnrpdb2_mibig_norine_path = (nerpa_dir
                                 / 'data'
                                 / 'input'
                                 / f'{pnrpdb_pref}_mibig_norine.tsv')
    is_mibig_norine_col = pnrpdb2_df['ID'].apply(is_mibig_norine)
    df_mibig_norine = pnrpdb2_df[is_mibig_norine_col]
    df_mibig_norine.to_csv(pnrpdb2_mibig_norine_path, sep='\t', index=False)

    # 2. pnrpdb2_deduplicated.tsv
    pnrpdb2_deduplicated_path = (nerpa_dir
                                 / 'data'
                                 / 'input'
                                 / f'{pnrpdb_pref}_deduplicated.tsv')
    is_representative_col = pnrpdb2_df['ID'].apply(is_representative)
    df_representatives = pnrpdb2_df[is_representative_col]
    df_representatives.to_csv(pnrpdb2_deduplicated_path, sep='\t', index=False)

    # 3. pnrpdb2_mibig_norine_deduplicated.tsv
    pnrpdb2_mibig_norine_deduplicated_path = (nerpa_dir
                                              / 'data'
                                              / 'input'
                                              / f'{pnrpdb_pref}_mibig_norine_deduplicated.tsv')
    df_mibig_norine_deduplicated = pnrpdb2_df[is_mibig_norine_col & is_representative_col]
    df_mibig_norine_deduplicated.to_csv(pnrpdb2_mibig_norine_deduplicated_path, sep='\t', index=False)

    # 4. Dump preprocessed NRP variants and rBAN records
    preprocessed_dir = nerpa_dir / 'data' / 'input' / 'preprocessed'
    preprocessed_dir.mkdir(exist_ok=True, parents=True)
    for table_name in (
            f'{pnrpdb_pref}',
            f'{pnrpdb_pref}_mibig_norine',
            f'{pnrpdb_pref}_deduplicated',
            f'{pnrpdb_pref}_mibig_norine_deduplicated',
    ):
        table = nerpa_dir / 'data' / 'input' / f'{table_name}.tsv'
        table_df = pd.read_csv(table, sep='\t')
        nrp_ids_in_table = set(table_df['ID'])

        rban_records_out_path = preprocessed_dir / f'{table_name}_parsed_rban_records.yaml'

        rban_records_subset = [
            rban_record.to_dict()
            for rban_record in rban_records
            if rban_record.compound_id in nrp_ids_in_table
        ]

        with open(rban_records_out_path, 'w') as f:
            yaml.dump(rban_records_subset, f)


if __name__ == '__main__':
    main()
