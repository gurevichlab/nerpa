import subprocess
import os
import sys
from pathlib import Path
import polars as pl
import yaml
from typing import Dict, List
import argparse
from joblib import Parallel, delayed

from src.benchmarking.data_frames import MIBiG_BGCs_Info


def parse_args(nerpa_dir: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="run Nerpa on training BGCs against their NRPs -- one Nerpa run per BGC")
    parser.add_argument(
        "-o", "--output-dir",
        dest="output_dir",
        type=Path,
        help=f"Output directory"
    )
    parser.add_argument(
        "-t",
        "--num-threads",
        dest="num_threads",
        type=int,
        default=1,
        help="Number of parallel jobs to run (passed to joblib.Parallel; -1 uses all cores)",
    )

    parser.add_argument(
        "--antismash-results",
        type=Path,
        help=(
            "A file with antismash paths for the BGCs -- one path per line."
        )
    )

    parser.add_argument(
        "--mibig-bgcs-info",
        type=Path,
    )
    parser.add_argument(
        "--nrps",
        type=Path,
        help="Path to the NRP database TSV",
    )
    return parser.parse_args()


def run_nerpa_for_bgc_and_its_nrps(
        nerpa_root: Path,
        nerpa_script: Path,
        bgc_id: str,
        nrp_ids: List[str],
        nrps_df: pl.DataFrame,
        antismash_results: Path,
        main_output_dir: Path,
        idx: int,
        total: int
) -> None:
        print(f'{idx}/{total}. Matching BGC {bgc_id} against {len(nrp_ids)} NRPs:\n{nrp_ids}')

        # Write the relevant NRPs to a temporary TSV file
        nrps_table: pl.DataFrame = (
            nrps_df.filter(pl.col('ID').is_in(nrp_ids))
        )
        smiles_tsv_path: Path = nerpa_root / 'tmp' / f'{bgc_id}_nrps.tsv'
        smiles_tsv_path.parent.mkdir(parents=True, exist_ok=True)
        nrps_table.write_csv(smiles_tsv_path, separator='\t')

        output_dir = main_output_dir / bgc_id

        # Construct the command
        command = [
            "python3", str(nerpa_script),
            "--antismash", str(antismash_results),
            "--smiles-tsv", smiles_tsv_path,
            "--col-id", "ID",
            "--output-dir", str(output_dir),
            "--force-output-dir",
            "--let-it-crash",
            "--dump-all-preprocessed",
        ]

        try:
            subprocess.run(
                command,
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        except subprocess.CalledProcessError:
            print(f"{bgc_id} failed.", file=sys.stderr)
        else:
            print(f"{bgc_id} success.")
    

def main():
    nerpa_root = Path(__file__).resolve().parent.parent.parent
    assert nerpa_root.name.startswith('nerpa'), f"Expected nerpa_dir to be the nerpa repo, got {nerpa_root}"
    nerpa_script = nerpa_root / "nerpa.py"

    args = parse_args(nerpa_root)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    antismash_results: list[Path] = [
        Path(line)
        for line in args.antismash_results.read_text().split("\n")
    ]
    bgc_to_as_path: dict[str, Path] = {
        as_path.name: as_path
        for as_path in antismash_results
    }
    bgc_ids: set[str] = set(bgc_to_as_path.keys())

    mibig_bgcs_info = MIBiG_BGCs_Info.from_csv(args.mibig_bgcs_info)
    bgc_to_nrps: Dict[str, List[str]] = {
        bgc_id: nrp_ids.split(",")
        for bgc_id, nrp_ids in mibig_bgcs_info.select(
            "bgc_id", "nrp_ids"
        ).iter_rows()
        if bgc_id in bgc_ids
    }
    nrps_df: pl.DataFrame = pl.read_csv(args.nrps, separator="\t")

    Parallel(n_jobs=args.num_threads)(
        delayed(run_nerpa_for_bgc_and_its_nrps)(
            nerpa_root=nerpa_root,
            nerpa_script=nerpa_script,
            bgc_id=bgc_id,
            nrp_ids=bgc_to_nrps[bgc_id],
            nrps_df=nrps_df,
            antismash_results=bgc_to_as_path[bgc_id],
            main_output_dir=args.output_dir,
            idx=i,
            total=len(bgc_to_nrps)
        )
        for i, bgc_id in enumerate(bgc_ids, start=1)
    )


if __name__ == '__main__':
    main()
