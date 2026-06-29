import os
import sys

import pandas as pd
import json
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple
import argparse

import yaml

from src.rban_parsing.rban_parser import Parsed_rBAN_Record


def run_nerpa_on_pnrpdb2(
        nerpa_dir: Path,
        pnrpdb2_path: Path,
        output_dir: Path,
        deduplicate: bool = True,
):
    output_dir.mkdir(parents=True, exist_ok=True)
    nerpa_script = nerpa_dir / "nerpa.py"

    num_threads = 32 if os.cpu_count() >= 40 else 8  # to work on cluster and locally

    # Any antismash results will do -- they are irrelevant but required by nerpa.py. 
    antismash_results = nerpa_dir / "test_data" / "antismash"

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
        "--dump-all-preprocessed",
        # "--keep-intermediate-files",
    ]

    if not deduplicate:
        command += ["--disable-deduplication"]

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
    nrp_representatives_path = (
        nerpa_results
        / 'preprocessed_input'
        / 'NRP_representatives.yaml'
    )
    with nrp_representatives_path.open() as f:
        nrp_representatives = yaml.safe_load(f)
    nrp_id_to_repr_id = {entry['nrp_id'][0]: entry['representative_id'][0]
                         for entry in nrp_representatives}
    return nrp_id_to_repr_id


def parse_args(nerpa_dir: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Preprocess NRP database by Nerpa.")

    default_nrp_database = nerpa_dir / 'data' / 'input' / 'pnrpdb2.tsv'
    parser.add_argument(
        "--nrp-database",
        type=Path,
        default=default_nrp_database,
        help=f"Path to a TSV file with NRPs to preprocess (if not provided, will use {default_nrp_database})."
    )

    default_preprocessed_dir = nerpa_dir / 'data' / 'input' / 'preprocessed'
    parser.add_argument(
        "--preprocessed-dir",
        type=Path,
        default=default_preprocessed_dir,
        help=f"Path to the directory where preprocessed data will be stored (default: {default_preprocessed_dir})."
    )

    default_filtered_tables_dir = nerpa_dir / 'data' / 'input'
    parser.add_argument(
        "--filtered-tables-dir",
        type=Path,
        default=default_filtered_tables_dir,
        help=f"Path to the directory where filtered tables will be stored (default: {default_filtered_tables_dir})."
    )

    modes = parser.add_subparsers(dest='mode', required=True)

    reuse_nerpa_results_parser = modes.add_parser(
        'reuse-nerpa-results',
        help="use existing Nerpa run on the database compouns"
    )
    reuse_nerpa_results_parser.add_argument(
        "--nerpa-results-dir",
        type=Path,
        required=True,
        help="Path to the directory containing the results of a previous Nerpa run on the database",
    )

    run_nerpa = modes.add_parser(
        'run-nerpa',
        help="run Nerpa on the database compounds"
    )

    default_nerpa_results_output_dir = (
        nerpa_dir
        / 'nerpa_results'
        / 'preprocess_pnrpdb2'
    )
    run_nerpa.add_argument(
        "--output-for-nerpa-run",
        type=Path,
        default=default_nerpa_results_output_dir,
        help=f"Path to the directory where Nerpa results will be stored (default: {default_nerpa_results_output_dir})."
    )

    run_nerpa.add_argument(
        "--deduplicate",
        action='store_true',
        help="Whether to run nerpa with deduplication enabled (default: True). If False, will disable deduplication in nerpa. Note that deduplication takes a long time, so you may want to disable it for testing purposes. However, the final results tables will be more useful with deduplication enabled."
    )

    args = parser.parse_args()
    return args


def main():
    nerpa_dir = Path(__file__).resolve().parent.parent
    assert (nerpa_dir / 'nerpa.py').exists(), f"nerpa.py not found in {nerpa_dir}"

    args = parse_args(nerpa_dir)
    pnrpdb2_path = args.nrp_database
    pnrpdb_pref = (
        pnrpdb2_path.stem[:-len('_raw')]
        if pnrpdb2_path.stem.endswith('_raw')
        else pnrpdb2_path.stem + '_filtered'
    )
    print(f'Using pnrpdb2 path: {pnrpdb2_path}, prefix for output tables: {pnrpdb_pref}')
    pnrpdb2_df = pd.read_csv(pnrpdb2_path, sep='\t')

    if args.mode == 'reuse-nerpa-results':
        nerpa_results = args.nerpa_results_dir
        print(f"Reusing existing Nerpa results from {nerpa_results}")
    else:
        nerpa_results = args.output_for_nerpa_run
        print(f"Running Nerpa on {pnrpdb2_path} and storing results in {args.output_for_nerpa_run}")
        
        run_nerpa_on_pnrpdb2(
            nerpa_dir=nerpa_dir,
            pnrpdb2_path=pnrpdb2_path,
            output_dir=args.output_for_nerpa_run,
            deduplicate=args.deduplicate,
        )

    nrp_id_to_repr_id = load_nrp_to_representative(nerpa_results)
    rban_records = [
        Parsed_rBAN_Record.from_dict(record)
        for record in yaml.safe_load((nerpa_results / 'preprocessed_input' / 'parsed_rban_records.yaml').open())
    ]
    processed_ids = set(nrp_id_to_repr_id.keys())  # IDs that Nerpa was able to process
    pnrpdb2_df = pnrpdb2_df[pnrpdb2_df['ID'].isin(processed_ids)]

    def is_mibig_norine(nrp_id: str) -> bool:
        return nrp_id.startswith('BGC') or nrp_id.startswith('NOR')

    def is_representative(nrp_id: str) -> bool:
        return nrp_id in nrp_id_to_repr_id.values()

    tables: List[Tuple[str, pd.DataFrame]] = []

    # 0. The filtered pnrpdb2 table with only records that were processed by nerpa
    tables.append((pnrpdb_pref, pnrpdb2_df))

    # 1. Only the subset of pnrpdb2 that are in MIBiG or Norine
    is_mibig_norine_col = pnrpdb2_df['ID'].apply(is_mibig_norine)
    df_mibig_norine = pnrpdb2_df[is_mibig_norine_col]
    tables.append((f'{pnrpdb_pref}_mibig_norine', df_mibig_norine))

    if args.deduplicate:
        # 2. Only the representatives from nerpa deduplication
        is_representative_col = pnrpdb2_df['ID'].apply(is_representative)
        df_representatives = pnrpdb2_df[is_representative_col]
        tables.append((f'{pnrpdb_pref}_deduplicated', df_representatives))

        # 3. Only the MIBiG or Norine NRPs with distinct representatives (i.e. deduplicated)
        repr_ids = df_mibig_norine['ID'].map(nrp_id_to_repr_id)
        df_mibig_norine_deduplicated = (
            df_mibig_norine
            .assign(_repr_id=repr_ids)
            .drop_duplicates(subset='_repr_id', keep='first')
            .drop(columns=['_repr_id'])
        )
        tables.append((f'{pnrpdb_pref}_mibig_norine_deduplicated', df_mibig_norine_deduplicated))

    # 4. Dump tables and preprocessed rBAN records
    for table_name, table_df in tables:
        nrp_ids_in_table = set(table_df['ID'])

        # Dump the table
        tsv_path = args.filtered_tables_dir / f'{table_name}.tsv'
        table_df.to_csv(tsv_path, sep='\t', index=False)

        # Dump the subset of rBAN records corresponding to the NRP IDs in this table
        rban_records_out_path = args.preprocessed_dir / f'{table_name}_parsed_rban_records.yaml'

        rban_records_subset = [
            rban_record.to_dict()
            for rban_record in rban_records
            if rban_record.compound_id in nrp_ids_in_table
        ]

        with open(rban_records_out_path, 'w') as f:
            yaml.dump(rban_records_subset, f)


if __name__ == '__main__':
    main()
