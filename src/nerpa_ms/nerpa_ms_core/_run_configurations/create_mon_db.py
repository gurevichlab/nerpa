# ARGS=""
import subprocess
import os
import argparse
from pathlib import Path

nerpa_ms_dir = Path(__file__).parent.parent

nerpa_dir = nerpa_ms_dir
while nerpa_dir.name != "src":
    nerpa_dir = nerpa_dir.parent
nerpa_dir = nerpa_dir.parent

parser = argparse.ArgumentParser(description="Run create_monomers_db with input and output paths")
parser.add_argument("--parsed-rban-records", "-i",
                    default=nerpa_dir / "nerpa_results/test_run/preprocessed_input/parsed_rban_records.yaml",
                    help="Path to YAML with the parsed rBAN records")
parser.add_argument("--out", "-o",
                    default=nerpa_ms_dir / "data/monomers_db.json",
                    help="Path to write the monomers DB JSON output")
args = parser.parse_args()

# Run the create_monomers_db binary with the provided paths
os.environ["RUST_BACKTRACE"] = "1"
subprocess.run([
    "cargo", "run", "--bin", "create_monomers_db", "--",
    "--parsed-rban-records", args.parsed_rban_records,
    "--out", args.out
])

