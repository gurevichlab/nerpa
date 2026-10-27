import os
import sys

import polars as pl
import json
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, NamedTuple, Optional, Callable, Hashable
import argparse

import yaml

from src.rban_parsing.rban_parser import Parsed_rBAN_Record
import re
import shutil


def run_nerpa_on_pnrpdb2(
        nerpa_dir: Path,
        pnrpdb2_path: Path,
        output_dir: Path,
        deduplicate: bool = True,
) -> None:
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

def load_nrp_to_representative(
        nerpa_results: Path,
) -> Dict[str, str]:
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


def rotate_represantative_mapping(
        nrp_id_to_repr_id: Dict[str, str],
        allowed_representatives: set[str],
) -> Dict[str, str]:
    """
    Given a mapping from NRP IDs to their representative IDs, and a set of allowed representatives,
    this function rotates the mapping so that each NRP ID maps to the first representative in its
    isomorphic group that is also in the allowed representatives set. If no representative in the
    group is allowed, an error is raised. The function returns a new mapping from NRP IDs to their allowed representative IDs.
    """
    # Create a reverse mapping from representative IDs to their corresponding NRP IDs
    repr_id_to_nrp_ids: Dict[str, List[str]] = {}
    for nrp_id, repr_id in nrp_id_to_repr_id.items():
        if repr_id not in repr_id_to_nrp_ids:
            repr_id_to_nrp_ids[repr_id] = []
        repr_id_to_nrp_ids[repr_id].append(nrp_id)

    # Create a new mapping from NRP IDs to their allowed representative IDs
    rotated_mapping: Dict[str, str] = {}
    for repr_id, nrp_ids in repr_id_to_nrp_ids.items():
        # Find the first allowed representative in this group
        allowed_repr = next((r for r in nrp_ids if r in allowed_representatives), None)
        if allowed_repr:
            for nrp_id in nrp_ids:
                rotated_mapping[nrp_id] = allowed_repr
        else:
            raise ValueError(f"No allowed representative found for the isomorphic group of representative ID '{repr_id}'.")
            

    return rotated_mapping


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


def expand_pnrpdb2_df(pnrpdb2_df: pl.DataFrame) -> pl.DataFrame:
    """
    PNRPDB2 rows have a column "Alternative IDs" which contains a comma-separated list of alternative IDs for the same compound.
    This function expands the DataFrame so that each alternative ID gets its own row, while keeping the other columns the same.
    """
    alternative_ids = (
        pl.col("Alternative IDs")
        .fill_null("")
        .cast(pl.String)
        .str.extract_all(r"[^,;| ]+")
    )

    return (
        pnrpdb2_df
        .with_columns(
            pl.concat_list(
                pl.col("ID").cast(pl.String),
                alternative_ids,
            ).alias("_expanded_id")
        )
        .explode("_expanded_id")
        .with_columns(pl.col("_expanded_id").alias("ID"))
        .drop("_expanded_id")
        .sort("ID")
    )


class TableData(NamedTuple):
    name: str
    df: pl.DataFrame
    readme: str


def handle_pnrpdb_expanded(
        pnrpdb2_raw_df: pl.DataFrame,
        tables_dir: Path,
        preprocessed_dir: Path,
        table_name_prefix: str,
        nerpa_output_dir: Path,
        nerpa_dir: Path,
) -> None:
    pnrpdb2_raw_expanded_df = expand_pnrpdb2_df(pnrpdb2_raw_df)
    pnrpdb2_raw_expanded_path = tables_dir / f'{table_name_prefix}_expanded.tsv'
    pnrpdb2_raw_expanded_df.write_csv(
        pnrpdb2_raw_expanded_path,
        separator='\t',
    )

    with open(pnrpdb2_raw_expanded_path.with_suffix('.tsv.readme'), 'w') as f:
        f.write(
            "This file is an expanded version of the original PNRPDB2 raw data, where each alternative ID has its own row.\n"
            "The 'ID' column contains the alternative ID, while other columns are copied from the original row.\n"
            "The original 'Alternative IDs' column is preserved in each row."
        )

    nerpa_output_dir_expanded = nerpa_output_dir.with_name(
        f"{nerpa_output_dir.name}_expanded"
    )
    print(
        f"Running Nerpa on {pnrpdb2_raw_expanded_path} "
        f"and storing results in {nerpa_output_dir_expanded}"
    )

    run_nerpa_on_pnrpdb2(
        nerpa_dir=nerpa_dir,
        pnrpdb2_path=pnrpdb2_raw_expanded_path,
        output_dir=nerpa_output_dir_expanded
    )

    # copy the parsed_rban_records.yaml from the nerpa output to the preprocessed directory
    source_path = (
        nerpa_output_dir_expanded
        / "preprocessed_input"
        / "parsed_rban_records.yaml"
    )
    destination_path = (
        preprocessed_dir
        / f"{table_name_prefix}_expanded_parsed_rban_records.yaml"
    )
    shutil.copy2(source_path, destination_path)


def main() -> None:
    nerpa_dir = Path(__file__).resolve().parent.parent.parent
    assert (nerpa_dir / 'nerpa.py').exists(), f"nerpa.py not found in {nerpa_dir}"

    args = parse_args(nerpa_dir)
    pnrpdb_pref = (
        args.pnrpdb2_raw.stem[:-len('_raw')]
        if args.pnrpdb2_raw.stem.endswith('_raw')
        else args.pnrpdb2_raw.stem + '_filtered'
    )
    print(f'Using pnrpdb2 path: {args.pnrpdb2_raw}, prefix for output tables: {pnrpdb_pref}')
    pnrpdb2_raw_df = pl.read_csv(args.pnrpdb2_raw, separator='\t')

    # pnrpdb2_expanded is different from other tables because it's not a filtered version of pnrpdb2_raw, but an expanded version of it.
    # So it was the least headache to handle it completely independently
    handle_pnrpdb_expanded(
        pnrpdb2_raw_df=pnrpdb2_raw_df,
        tables_dir=args.filtered_tables_dir,
        preprocessed_dir=args.preprocessed_dir,
        table_name_prefix=pnrpdb_pref,
        nerpa_output_dir=args.output_for_nerpa_run,
        nerpa_dir=nerpa_dir,
    )

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

    # pnrpdb2_raw_nrps = set(pnrpdb2_raw_df.get_column("ID").to_list())
    nrp_id_to_repr_id: Dict[str, str] = load_nrp_to_representative(nerpa_results)
    # nrp_id_to_repr_id = rotate_represantative_mapping(nrp_id_to_repr_id, pnrpdb2_raw_nrps)

    rban_records: List[Parsed_rBAN_Record] = [
        Parsed_rBAN_Record.from_dict(record)
        for record in yaml.safe_load(
            (nerpa_results / "preprocessed_input" / "parsed_rban_records.yaml").open()
        )
    ]
    processed_ids: set[str] = set(nrp_id_to_repr_id)
    pnrpdb2_df: pl.DataFrame = pnrpdb2_raw_df.filter(
        pl.col("ID").is_in(list(processed_ids))
    )

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
    df_mibig_norine = pnrpdb2_df.filter(
        pl.col("ID").str.starts_with("BGC")
        | pl.col("ID").str.starts_with("NOR")
    )
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
    representative_ids = list(nrp_id_to_repr_id.values())
    df_representatives = pnrpdb2_df.filter(
        pl.col("ID").is_in(representative_ids)
    )
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
    df_mibig_norine_deduplicated = (
        df_mibig_norine
        .with_columns(
            pl.col("ID")
            .replace_strict(nrp_id_to_repr_id)
            .alias("_repr_id")
        )
        .unique(subset="_repr_id", keep="first", maintain_order=True)
        .drop("_repr_id")
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
        nrp_ids_in_table = set(table.df.get_column("ID").to_list())

        # Dump the table
        tsv_path = args.filtered_tables_dir / f'{table.name}.tsv'
        table.df.write_csv(tsv_path, separator='\t')

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
