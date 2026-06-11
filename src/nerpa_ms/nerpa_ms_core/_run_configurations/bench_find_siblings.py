# ARGS=""
import subprocess
import os
from pathlib import Path
# os.environ["RUST_BACKTRACE"] = "1"
os.environ["RUSTFLAGS"] = "-Awarnings"


nerpa_root = (
    Path(__file__).resolve()
    .parent.parent.parent.parent.parent
)
assert (nerpa_root / "nerpa.py").exists(), f"Invalid nerpa_root: {nerpa_root}"

nerpa_results = nerpa_root / "nerpa_results/mibig_bgc_vs_its_nrps"
compound_distances = nerpa_root / "data/for_training_and_testing/nrp_distances.txt"

subprocess.run(["cargo", "build", "-q"])
command = ' '.join(['cargo run --bin benchmark_find_siblings --',
                    f'--nerpa-results {nerpa_results}',
                    # f'--compound-distances {compound_distances}',
                    '--max-edits 1',
                    '--num-variants-per-num-edits 1000',
                    '--out benchmarking/find_siblings',
                    '--nerpa-root', str(nerpa_root),
                    '--monomers-db-json data/monomers_db.json',
                    ])
subprocess.run(command.split())

