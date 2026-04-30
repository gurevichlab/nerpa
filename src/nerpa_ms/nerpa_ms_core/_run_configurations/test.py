# ARGS=""
import subprocess
import os
from pathlib import Path
os.environ["RUST_BACKTRACE"] = "0"
os.environ["RUSTFLAGS"] = "-Awarnings"


nerpa_root = (
    Path(__file__).resolve()
    .parent.parent.parent.parent.parent
)
assert (nerpa_root / "nerpa.py").exists(), f"Invalid nerpa_root: {nerpa_root}"

subprocess.run(["cargo", "build", "-q"])
command = ' '.join(['cargo run --bin nerpa_ms_core --',
                    '--nerpa-results fixtures/nerpa_results_1214.5',
                    '--max-edits 1',
                    '--num-variants-per-num-edits 5',
                    '--out output',
                    '--monomers-db-json data/monomers_db.json',
                    '--nerpa-root', str(nerpa_root),])
subprocess.run(command.split())

