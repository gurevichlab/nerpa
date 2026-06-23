# ARGS=""
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
        "benchmark Nerpa-MS results on MIBiG BGCs against their NRPs. "
        "See how often Nerpa-MS generates an NRP corresponding to a BGC given its sibling as the template"
    )
    default_nerpa_results = (
        nerpa_root
        / 'nerpa_results'
        / 'all_mibig_bgcs_vs_their_nrps'
    )
    parser.add_argument(
        "--nerpa-results",
        dest="nerpa_results",
        default=default_nerpa_results,
        type=Path,
        help=(
            "Path to the nerpa_results directory. "
            "Each subdirectory of this directory should be a Nerpa run "
            "on a MIBiG BGC vs its NRP(s), named after the BGC ID (e.g. BGC0001234). "
            f"Default: {default_nerpa_results}"
        ),
    )
    default_nerpa_ms_results = (
        nerpa_ms_root
        / 'output'
        / 'benchmarking'
        / 'find_siblings'
        / 'nerpa_ms_results'
    )
    parser.add_argument(
        "--nerpa-ms-results",
        dest="nerpa_ms_results",
        default=default_nerpa_ms_results,
        type=Path,
        help=(
            "Path to the nerpa__ms_results directory. "
            "Each subdirectory of this directory should be a Nerpa-MS run "
            "on a MIBiG BGC vs its NRP(s), named after the BGC ID (e.g. BGC0001234). "
            f"Default: {default_nerpa_ms_results}"
        ),
    )

    default_out_dir = (
        nerpa_ms_root
        / 'output'
        / 'benchmarking'
        / 'find_siblings'
    )
    parser.add_argument(
        "-o", "--output-dir",
        dest="output_dir",
        default=default_out_dir,
        type=Path,
        help=f"Output directory (default: {default_out_dir}"
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
        'cargo', 'run', '--bin', 'benchmark_find_siblings', '--',
        f'--nerpa-results={args.nerpa_results}',
        f'--nerpa-ms-results={args.nerpa_ms_results}',
        f'--out={args.output_dir}',
        '--nerpa-root', str(nerpa_root),
        '--monomers-db-json', 'data/monomers_db.json',
    ]
    print(f"Running in {nerpa_ms_root}:\n  {' '.join(command)}")
    subprocess.run(command, cwd=nerpa_ms_root, check=True)

if __name__ == "__main__":
    main()
