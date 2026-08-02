import subprocess
import os
from pathlib import Path
import argparse

def parse_args(
        nerpa_root: Path,
        nerpa_ms_root: Path
):
    parser = argparse.ArgumentParser(
        prog=
        "plot benchmarking of Nerpa-MS results on MIBiG BGCs against their NRPs. "
        "See how often Nerpa-MS generates an NRP corresponding to a BGC given its sibling as the template"
    )
    default_benchmark_results = (
        nerpa_ms_root
        / 'output'
        / 'benchmarking'
        / 'find_siblings'
        / 'benchmarking_results.json'
    )
    parser.add_argument(
        "--benchmarking-results",
        dest="benchmarking_results",
        default=default_benchmark_results,
        type=Path,
        help=(
            "Path to a JSON with the benchmarking results (find siblings test). "
            f"Default: {default_benchmark_results}"
        ),
    )
    default_out_file = (
        nerpa_ms_root
        / 'output'
        / 'benchmarking'
        / 'find_siblings'
        / 'plot.svg'
    )
    parser.add_argument(
        "-o", "--out",
        dest="out",
        default=default_out_file,
        type=Path,
        help=f"Output file (default: {default_out_file}"
    )

    return parser.parse_args()

def main():
    os.environ["RUST_BACKTRACE"] = "1"
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
    args = parse_args(nerpa_root=nerpa_root, nerpa_ms_root=nerpa_ms_root)
    command = [
        'cargo', 'run', '--bin', 'benchmark_find_siblings_plot', '--',
        f'--target-search-results={args.benchmarking_results}',
        f'--out={args.out}',
        '--nerpa-root', str(nerpa_root),
    ]
    print(f"Running in {nerpa_ms_root}:\n  {' '.join(command)}")
    subprocess.run(command, cwd=nerpa_ms_root, check=True)

if __name__ == "__main__":
    main()
