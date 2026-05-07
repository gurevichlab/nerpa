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
                    '--nerpa-results fixtures/test_run',
                    '--max-nerpa-matches 3',
                    '--max-edits 2',
                    '--num-variants-per-num-edits 3',
                    '--out output',
                    '--monomers-db-json data/monomers_db.json',
                    '--nerpa-root', str(nerpa_root),])
subprocess.run(command.split())

