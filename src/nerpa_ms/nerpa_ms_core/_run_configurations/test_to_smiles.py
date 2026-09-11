# ARGS=""
import subprocess
import os
from pathlib import Path

os.environ["RUST_BACKTRACE"] = "1"

nerpa_ms_root = Path(__file__).parent.parent
nerpa_root = nerpa_ms_root.parent.parent.parent

if not (nerpa_root / "nerpa.py").exists():
    raise FileNotFoundError(f"Invalid nerpa_root: {nerpa_root}")

rban_records_path = (
    nerpa_root
    / "data/input/preprocessed"
    / "pnrpdb2_mibig_norine_parsed_rban_records.yaml"
)

command = ' '.join([
    'cargo run --bin test_to_smiles --',
    f'--parsed-rban-records {rban_records_path}',
])
subprocess.run(command.split())

