import subprocess
import os
from pathlib import Path
import argparse

def get_default_out_dir(
        nerpa_ms_root: Path,
        input_dir_name: str
) -> Path:
    return (
        nerpa_ms_root
        / 'output'
        / input_dir_name
    )

def get_nerpa_results_dir(
        nerpa_resuls_all: Path,
        dir_substr: str
) -> Path:
    matches = [
        subdir
        for subdir in nerpa_resuls_all.iterdir()
        if dir_substr in subdir.name
    ]
    if len(matches) == 0:
        raise ValueError(f"No subdirectory of {nerpa_resuls_all} contains '{dir_substr}'")
    elif len(matches) > 1:
        raise ValueError(f"Multiple subdirectories of {nerpa_resuls_all} contain '{dir_substr}': {[match.name for match in matches]}")
    return matches[0]
    

def parse_args(
        nerpa_root: Path,
        nerpa_ms_root: Path
):
    parser = argparse.ArgumentParser(
        prog=
        "run Nerpa on MIBiG BGC with >= min_a_domains A domains "
        "vs its NRP(s) with >= min_monomers NRPS monomers"
    )
    default_max_edits = 3
    parser.add_argument(
        "--max-edits",
        dest="max_edits",
        type=int,
        default=default_max_edits,
        help=f"Maximum number of edits to consider when generating variants (default: {default_max_edits})",
    )

    default_num_variants_per_num_edits = 1000
    parser.add_argument(
        "--num-variants-per-num-edits",
        dest="num_variants_per_num_edits",
        type=int,
        default=default_num_variants_per_num_edits,
        help=f"Number of variants to generate per number of edits (default: {default_num_variants_per_num_edits})",
    )

    default_out_dir =  get_default_out_dir(
        nerpa_ms_root=nerpa_ms_root,
        input_dir_name='[INPUT_DIR_NAME]'
    )
        
    parser.add_argument(
        "-o", "--output-dir",
        dest="output_dir",
        type=Path,
        help=f"Output directory (default: {default_out_dir}"
    )
    default_nerpa_results = (
        nerpa_root
        / 'nerpa_results'
        / 'all_mibig_bgcs_vs_their_nrps'
    )
    parser.add_argument(
        "--nerpa-results-all",
        dest="nerpa_results_all",
        default=default_nerpa_results,
        type=Path,
        help=(
            "Path to the nerpa_results directory. "
            "Each subdirectory of this directory should be a Nerpa run "
            "on a MIBiG BGC vs its NRP(s), named after the BGC ID (e.g. BGC0001234). "
            f"Default: {default_nerpa_results}"
        ),
    )
    parser.add_argument(
        "--bgc-id",
        dest="bgc_id",
        help=(
            "The name of the subdirectory in nerpa-results to run nerpa_ms on. "
            "Unique substring is acceptable (e.g. 1234 instead of BGC0001234)"
        ),
    )
    args = parser.parse_args()
    args.nerpa_results = get_nerpa_results_dir(
        args.nerpa_results_all,
        dir_substr=args.bgc_id
    )
    if args.output_dir is None:
        args.output_dir = get_default_out_dir(
            nerpa_ms_root=nerpa_ms_root,
            input_dir_name=args.nerpa_results.name
        )

    return args
            

def main():
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
    os.environ["RUST_BACKTRACE"] = "0"
    os.environ["RUSTFLAGS"] = "-Awarnings"

    command = [
        'cargo',  'run',
        '--bin', 'nerpa_ms_core', '--',
        '--nerpa-results', f'{args.nerpa_results}',
        '--max-edits', f'{args.max_edits}',
        '--num-variants-per-num-edits', f'{args.num_variants_per_num_edits}',
        '--out', f'{args.output_dir}',
        '--monomers-db-json', 'data/monomers_db.json',
        '--nerpa-root', str(nerpa_root),
    ]
    print(f"Running in {nerpa_ms_root}:\n  {' '.join(command)}")
    subprocess.run(command, cwd=nerpa_ms_root, check=True)


if __name__ == '__main__':
    main()
