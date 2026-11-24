from pathlib import Path

import polars as pl
import pandas as pd
from collections import Counter
import yaml
from src.monomer_names_helper.monomer_names_helper import MonomerNamesHelper

import argparse

def parse_cli_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--norine-graphs",
        type=Path,
        required=True,
        help="Path to norine_graphs.csv",
    )
    parser.add_argument(
        "--monomer-names-table",
        type=Path,
        required=True,
        help="Path to monomer_names_table.tsv",
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Path to output TSV with monomer counts",
    )
    args = parser.parse_args()
    return args
    

def main() -> None:
    args = parse_cli_args()

    # Load the Norine monomers table
    norine_graphs_table = pd.read_csv(args.norine_graphs, delimiter=";")
    graph_records = norine_graphs_table["graph"].to_list()

    monomers_str = []
    for graph_record in graph_records:
        monomers_str.extend(graph_record.split("@")[0].split(","))

    # Count monomer cnts
    monomers_counter = Counter(monomers_str)

    monomer_names_table = pl.read_csv(args.monomer_names_table, separator="\t")
    monomer_names_helper = MonomerNamesHelper(names_table=monomer_names_table)

    rows = []
    for monomer_name, count in monomers_counter.items():
        monomer = monomer_names_helper.parsed_name(monomer_name, name_format="rBAN/Norine")
        rows.append({
            "Norine_name": monomer_name,
            "count": count,
            "Nerpa_residue": monomer.residue,
            "methylated": monomer.methylated,
            "chirality": monomer.chirality.name,
            "is_pks_hybrid": monomer.is_pks_hybrid,
        })

    # Write the output TSV
    out_df = pl.DataFrame(rows).sort("count", descending=True)
    out_df.write_csv(args.out, separator="\t")


if __name__ == "__main__":
    main()
