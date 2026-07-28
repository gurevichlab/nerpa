import os
import sys

import pandas as pd
import json
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, NamedTuple, Optional, Callable, Hashable
import argparse

import yaml

from src.rban_parsing.rban_parser import Parsed_rBAN_Record
import re


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

    default_nrp_database = nerpa_dir / 'data' / 'input' / 'pnrpdb2_raw.tsv'
    parser.add_argument(
        "--pnrpdb2-raw",
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

    default_nerpa_results_output_dir = (
        nerpa_dir
        / 'nerpa_results'
        / 'preprocess_pnrpdb2'
    )
    parser.add_argument(
        "--output-for-nerpa-run",
        type=Path,
        default=default_nerpa_results_output_dir,
        help=f"Path to the directory where Nerpa results will be stored (default: {default_nerpa_results_output_dir})."
    )

    args = parser.parse_args()
    return args


def expand_pnrpdb2_df(pnrpdb2_df: pd.DataFrame) -> pd.DataFrame:
    """
    PNRPDB2 rows have a column "Alternative IDs" which contains a comma-separated list of alternative IDs for the same compound.
    This function expands the DataFrame so that each alternative ID gets its own row, while keeping the other columns the same.
    """
    expanded_rows = []
    for _, row in pnrpdb2_df.iterrows():
        alt_ids = row.get('Alternative IDs', '')
        alt_ids_list = [
            s.strip()
            for s in re.split(r'[,;| ]+', alt_ids)  # treat any of "," ";", "|", " " as separators
            if s.strip()
        ]
        if not alt_ids_list:
            expanded_rows.append(row)
        else:
            for alt_id in alt_ids_list:
                new_row = row.copy()
                new_row['ID'] = alt_id
                expanded_rows.append(new_row)
    expanded_df = pd.DataFrame(expanded_rows)
    return expanded_df


class TableData(NamedTuple):
    name: str
    df: pd.DataFrame
    readme: str

def main():
    nerpa_dir = Path(__file__).resolve().parent.parent
    assert (nerpa_dir / 'nerpa.py').exists(), f"nerpa.py not found in {nerpa_dir}"

    args = parse_args(nerpa_dir)
    pnrpdb_pref = (
        args.pnrpdb2_raw.stem[:-len('_raw')]
        if args.pnrpdb2_raw.stem.endswith('_raw')
        else args.pnrpdb2_raw.stem + '_filtered'
    )
    print(f'Using pnrpdb2 path: {args.pnrpdb2_raw}, prefix for output tables: {pnrpdb_pref}')
    pnrpdb2_raw_df = pd.read_csv(args.pnrpdb2_raw, sep='\t')

    # pnrpdb2_raw_expanded_df = expand_pnrpdb2_df(pnrpdb2_raw_df)
    # pnrpdb2_raw_expanded_path = args.preprocessed_dir / f'{args.pnrpdb_raw.stem}_expanded.tsv'
    # with pnrpdb2_raw_expanded_path.open('w') as f:
    #     pnrpdb2_raw_expanded_df.to_csv(f, sep='\t', index=False)
    # with open(args.preprocessed_dir / f'{args.pnrpdb_raw.stem}_expanded.tsv.readme', 'w') as f:
    #     f.write(
    #         "This file is an expanded version of the original PNRPDB2 raw data, where each alternative ID has its own row.\n"
    #         "The 'ID' column contains the alternative ID, while other columns are copied from the original row.\n"
    #         "The original 'Alternative IDs' column is preserved in each row."
    #     )


    nerpa_results = args.output_for_nerpa_run
    print(
        f"Running Nerpa on {args.pnrpdb2_raw} "
        f"and storing results in {args.output_for_nerpa_run}"
    )

    run_nerpa_on_pnrpdb2(
        nerpa_dir=nerpa_dir,
        pnrpdb2_path=args.pnrpdb2_raw,
        output_dir=args.output_for_nerpa_run,
    )

    nrp_id_to_repr_id = load_nrp_to_representative(nerpa_results)
    rban_records = [
        Parsed_rBAN_Record.from_dict(record)
        for record in yaml.safe_load((nerpa_results / 'preprocessed_input' / 'parsed_rban_records.yaml').open())
    ]
    processed_ids = set(nrp_id_to_repr_id.keys())  # IDs that Nerpa was able to process
    pnrpdb2_df = pnrpdb2_raw_df[pnrpdb2_raw_df['ID'].isin(processed_ids)]

    def is_mibig_norine(nrp_id: str) -> bool:
        return nrp_id.startswith('BGC') or nrp_id.startswith('NOR')

    def is_representative(nrp_id: str) -> bool:
        return nrp_id in nrp_id_to_repr_id.values()

    tables: List[TableData] = []

    # 1. The filtered pnrpdb2 table with only records that were processed by nerpa
    tables.append(
        TableData(
            name=pnrpdb_pref,
            df=pnrpdb2_df,
            readme=(
                "This table contains only the records from PNRPDB2 for which Nerpa "
                "was able to generate NRP variants."
            )
        )
    )

    # 2. Only the subset of pnrpdb2 that are in MIBiG or Norine
    is_mibig_norine_col = pnrpdb2_df['ID'].apply(is_mibig_norine)
    df_mibig_norine = pnrpdb2_df[is_mibig_norine_col]
    tables.append(
        TableData(
            name=f'{pnrpdb_pref}_mibig_norine',
            df=df_mibig_norine,
            readme=(
                "This table contains only the records from PNRPDB2 "
                "which are \"confirmed NRPs\", i.e. originated from MIBiG or Norine."
            )
        )
    )

    # 3. Only the representatives from nerpa deduplication
    is_representative_col = pnrpdb2_df['ID'].apply(is_representative)
    df_representatives = pnrpdb2_df[is_representative_col]
    tables.append(
        TableData(
            name=f'{pnrpdb_pref}_deduplicated',
            df=df_representatives,
            readme=(
                "This table contains only the records from PNRPDB2 "
                "whose NRP Variants are representatives of their isomorphic groups."
            )
        )
    )

    # 3. Only the MIBiG or Norine NRPs with distinct representatives (i.e. deduplicated)
    repr_ids = df_mibig_norine['ID'].map(nrp_id_to_repr_id)
    df_mibig_norine_deduplicated = (
        df_mibig_norine
        .assign(_repr_id=repr_ids)
        .drop_duplicates(subset='_repr_id', keep='first')
        .drop(columns=['_repr_id'])
    )
    tables.append(
        TableData(
            name=f'{pnrpdb_pref}_mibig_norine_deduplicated',
            df=df_mibig_norine_deduplicated,
            readme=(
                "This table contains only the records from PNRPDB2 "
                "which are \"confirmed NRPs\" (from MIBiG or Norine) "
                "and whose NRP variants are representatives of their isomorphic groups."
            )
        )
    )

    # 4. Dump tables and preprocessed rBAN records
    for table in tables:
        nrp_ids_in_table = set(table.df['ID'])

        # Dump the table
        tsv_path = args.filtered_tables_dir / f'{table.name}.tsv'
        table.df.to_csv(tsv_path, sep='\t', index=False)

        # Dump the README for the table
        readme_path = args.filtered_tables_dir / f'{table.name}.tsv.readme'
        with open(readme_path, 'w') as f:
            f.write(table.readme)

        # Dump the subset of rBAN records corresponding to the NRP IDs in this table
        rban_records_out_path = args.preprocessed_dir / f'{table.name}_parsed_rban_records.yaml'

        rban_records_subset = [
            rban_record.to_dict()
            for rban_record in rban_records
            if rban_record.compound_id in nrp_ids_in_table
        ]

        with open(rban_records_out_path, 'w') as f:
            yaml.dump(rban_records_subset, f)

    # 5. Dump the unfiltered preprocessed rBAN records
    # (including the records that Nerpa was unable to process)
    all_rban_records_out_path = (
        args.preprocessed_dir
        / f'{args.pnrpdb2_raw.stem}_parsed_rban_records.yaml'
    )
    with open(all_rban_records_out_path, 'w') as f:
        yaml.dump([rban_record.to_dict() for rban_record in rban_records], f)


if __name__ == '__main__':
    main()
