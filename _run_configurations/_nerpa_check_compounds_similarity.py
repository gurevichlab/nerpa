#!/usr/bin/env python3
import argparse
import sys
from itertools import combinations
from pathlib import Path
import subprocess
import platform
import polars as pl
import yaml
from networkx.algorithms.isomorphism.isomorph import is_isomorphic

from src.generic.graphs import graphs_one_substitution_away
from src.rban_parsing.rban_parser import Parsed_rBAN_Record
from src.rban_parsing.nrp_variant_types import NRP_Variant
from scripts.pnrpdb_compound_similarity import nerpa_mon_cmp, unknown_chr_equal_known_cmp, parsed_record_to_graph

def usage() -> None:
    exe = Path(sys.argv[0]).name
    print(
        f"Usage:\n"
        f"  {exe} <id1> <id2> ...\n"
        f"  {exe} --file ids.txt\n"
        f"\n"
        f"Notes:\n"
        f"  - IDs must match the 'ID' column in ./data/input/pnrpdb2.tsv\n"
        f"  - Writes a filtered TSV and then runs nerpa for each ID.\n"
    )

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Read IDs from command line, file, or load-results directory."
    )

    parser.add_argument(
        "--ids",
        nargs="+",
        help="IDs provided directly",
    )

    parser.add_argument(
        "--file",
        type=Path,
        help="Path to file containing IDs (one per line, # for comments)",
    )

    parser.add_argument(
        "--load-results",
        help="Directory with parsed compound results",
    )

    return parser


def read_ids_from_file(p: Path) -> list[str]:
    if not p.exists():
        raise FileNotFoundError(f"IDs file not found: {p}")

    ids: list[str] = []
    for line in p.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        ids.append(line.split()[0])

    if not ids:
        raise ValueError(f"No IDs found in file: {p}")

    return ids


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = build_parser()
    args = parser.parse_args(argv[1:])
    if args.load_results == 'default':
        args.load_results = './nerpa_results/parsed_compounds/nerpa_results'
    return args


def open_dir(path: Path) -> None:
    path = path.resolve()
    try:
        subprocess.Popen(
            ["xdg-open", str(path)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            start_new_session=True
        )
    except Exception as e:
        print(f"Could not open directory automatically: {e}")
        print(f"Results are in: {path}")


def get_nerpa_results(args: argparse.Namespace,
                      nerpa_dir: Path,
                      pnrpdb_path: Path,
                      antismash_dir: Path) -> Path:
    if args.file:
        ids = read_ids_from_file(args.file)
    else:
        ids = [a.strip("',()") for a in args.ids if a.strip()]
        if not ids:
            raise ValueError("No IDs provided.")

    # Output layout
    output_root = nerpa_dir / "nerpa_results" / "parsed_compounds"
    output_root.mkdir(parents=True, exist_ok=True)

    # Read and filter DB
    pnrpdb = pl.read_csv(pnrpdb_path, separator="\t")
    missing_cols = [c for c in ("ID", "SMILES") if c not in pnrpdb.columns]
    if missing_cols:
        raise ValueError(f"pnrpdb2.tsv missing required columns: {missing_cols}")

    ids_set = set(ids)
    filtered = pnrpdb.filter(pl.col("ID").is_in(list(ids_set)))

    # Warn about IDs not found
    found_ids = set(filtered.select("ID").unique().to_series().to_list())
    not_found = [i for i in ids if i not in found_ids]
    if not_found:
        print("WARNING: Some IDs were not found in pnrpdb2.tsv:")
        for i in not_found:
            print(f"  - {i}")

    if filtered.height == 0:
        raise ValueError("None of the provided IDs were found; nothing to do.")

    # Write filtered TSV
    filtered_tsv = output_root / "compounds_filtered.tsv"
    filtered_tsv.parent.mkdir(parents=True, exist_ok=True)
    filtered.write_csv(filtered_tsv, separator="\t")
    print(f"Wrote filtered compounds TSV: {filtered_tsv}")


    nerpa_results_dir = output_root / 'nerpa_results'
    nerpa_script = nerpa_dir / "nerpa.py"

    # Run nerpa per compound
    command = [
        "python3",
        str(nerpa_script),
        "--antismash",
        str(antismash_dir),
        "--smiles-tsv",
        str(filtered_tsv),
        "--col-id", "ID",
        "--output-dir",
        str(nerpa_results_dir),
        "--force-output-dir",
        "--fast-matching",
        "--let-it-crash",
        "--dump-all-preprocessed",
    ]

    subprocess.run(command, check=True)
    return nerpa_results_dir


def main() -> None:
    nerpa_dir = Path(__file__).resolve().parent

    pnrpdb_path = nerpa_dir / "data/input/pnrpdb2.tsv"
    antismash_dir = nerpa_dir / "test_data/antismash"

    if not pnrpdb_path.exists():
        raise FileNotFoundError(f"pnrpdb2.tsv not found at: {pnrpdb_path}")
    if not antismash_dir.exists():
        raise FileNotFoundError(f"test antismash dir not found at: {antismash_dir}")

    args = parse_args(sys.argv)
    if args.load_results:
        nerpa_results_dir = Path(args.load_results).resolve()
    else:
        nerpa_results_dir = get_nerpa_results(args,
                                              nerpa_dir=nerpa_dir,
                                              pnrpdb_path=pnrpdb_path,
                                              antismash_dir=antismash_dir)


    # # Load parsed_rban_records.yaml
    # with open(nerpa_results_dir / 'preprocessed_input' / 'parsed_rban_records.yaml') as f:
    #     parsed_rban_records = [Parsed_rBAN_Record.from_dict(record)
    #                            for record in yaml.safe_load(f)]

    # Load NRP variants
    with open(nerpa_results_dir / 'preprocessed_input' / 'NRP_variants.yaml') as f:
        nrp_variants = [NRP_Variant.from_yaml_dict(record)
                        for record in yaml.safe_load(f)]

    # rban_records_by_id = {r.compound_id: r for r in parsed_rban_records}
    nrp_variants_by_id = {v.nrp_variant_id.nrp_id: v for v in nrp_variants}

    for id1, id2 in combinations(nrp_variants_by_id, 2):
        # rban_graph1 = parsed_record_to_graph(rban_records_by_id[id1])
        # rban_graph2 = rban_records_by_id[id2]
        nerpa_graph1 = nrp_variants_by_id[id1].to_nx_digraph()
        nerpa_graph2 = nrp_variants_by_id[id2].to_nx_digraph()
        are_isomorphic = is_isomorphic(nerpa_graph1, nerpa_graph2,
                                       node_match=lambda n1, n2: unknown_chr_equal_known_cmp(n1['monomer'], n2['monomer']))

        one_sub_away = graphs_one_substitution_away(
            nerpa_graph1, nerpa_graph2,
            nodes_comparator=unknown_chr_equal_known_cmp,
            label_key="monomer"
        )

        if are_isomorphic:
            print(f'{id1} and {id2} are isomorphic!')
        elif one_sub_away:
            print(f'{id1} and {id2} are one substitution away!')
        else:
            print(f'{id1} and {id2} are different.')


if __name__ == "__main__":
    main()
