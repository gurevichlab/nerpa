#!/usr/bin/env python3
import sys
from pathlib import Path
import subprocess
import platform
import polars as pl


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


def read_ids(argv: list[str]) -> list[str]:
    if len(argv) < 2:
        usage()
        raise SystemExit(2)

    if argv[1] in ("-h", "--help"):
        usage()
        raise SystemExit(0)

    if argv[1] == "--file":
        if len(argv) != 3:
            usage()
            raise SystemExit(2)
        p = Path(argv[2]).expanduser()
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

    # Otherwise treat remaining args as IDs
    ids = [a.strip("',()") for a in argv[1:] if a.strip()]
    if not ids:
        raise ValueError("No IDs provided.")
    return ids


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


def main() -> None:
    nerpa_dir = Path(__file__).resolve().parent.parent
    nerpa_script = nerpa_dir / "nerpa.py"

    pnrpdb_path = nerpa_dir / "data/input/pnrpdb2.tsv"
    antismash_dir = nerpa_dir / "test_data/antismash"

    if not nerpa_script.exists():
        raise FileNotFoundError(f"nerpa.py not found at: {nerpa_script}")
    if not pnrpdb_path.exists():
        raise FileNotFoundError(f"pnrpdb2.tsv not found at: {pnrpdb_path}")
    if not antismash_dir.exists():
        raise FileNotFoundError(f"test antismash dir not found at: {antismash_dir}")

    ids = read_ids(sys.argv)

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

    # Summary + open results folder
    print("\n=== Done ===")
    print(f"Results root: {output_root}")
    open_dir(nerpa_results_dir / 'NRP_images')


if __name__ == "__main__":
    main()
