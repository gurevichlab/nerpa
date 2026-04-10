# ARGS=""
import subprocess
import os
import argparse
from pathlib import Path

nerpa_ms_dir = Path(__file__).parent.parent

parser = argparse.ArgumentParser(description="Generate test inputs for dynamic programming")
parser.add_argument("--out", "-o",
                    default=nerpa_ms_dir / "tests/fixtures/dp_vs_bruteforce_tests.json",
                    help="Path to write the tests JSON output")
args = parser.parse_args()

os.environ["RUST_BACKTRACE"] = "1"
subprocess.run(["cargo", "build"])
subprocess.run([
    "cargo", "run", "--bin", "gen_dp_tests", "--",
    "--out", args.out,
    "--n", "50",
    "--seed", "1",
    "--max-weight", "2",
    "--alphabet",  "4"
])

