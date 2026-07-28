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
nerpa_ms_root = (
    nerpa_root
    / 'src'
    / 'nerpa_ms'
    / 'nerpa_ms_core'
)

subprocess.run(["cargo", "build"], cwd=nerpa_ms_root, check=True)
command = [
    'cargo',  'run',
    '--bin', 'nerpa_ms_core', '--',
    '--nerpa-results', '/home/ilianolhin/work/tools/nerpa/nerpa_results/myxochromides_wo_A',
    '--max-edits', '2',
    '--num-variants-per-num-edits', '1000',
    '--out', 'myxochromides_wo_A',
    '--monomers-db-json', 'data/monomers_db.json',
    '--nerpa-root', str(nerpa_root),
]
print(f"Running in {nerpa_ms_root}:\n  {' '.join(command)}")
subprocess.run(command, cwd=nerpa_ms_root, check=True)

