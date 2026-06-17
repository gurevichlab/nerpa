import subprocess
import os
import sys
from pathlib import Path
import polars as pl
import yaml
from typing import Dict, List
from src.benchmarking.data_frames import (
    MIBiG_BGCs_Info,
    PNRPDB_Info,
    get_good_mibig_bgs_to_nrps,
)
import argparse
from joblib import Parallel, delayed


def parse_args(nerpa_dir: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="run Nerpa on MIBiG BGC with >= min_a_domains A domains vs its NRP(s) with >= min_monomers NRPS monomers")
    default_output_dir = nerpa_dir / "nerpa_results/mibig_good_bgcs_vs_their_good_nrps"
    parser.add_argument(
        "-o", "--output-dir",
        dest="output_dir",
        type=Path,
        default=default_output_dir,
        help=f"Output directory (default: {default_output_dir}"
    )
    parser.add_argument(
        "--min-a-domains",
        dest="min_a_domains",
        type=int,
        default=3,
        help="Minimum number of A domains in the BGC to include it (default: 3)")

    parser.add_argument(
        "--min-monomers",
        dest="min_monomers",
        type=int,
        default=3,
        help="Minimum number of monomers in the NRP to include it (default: 3)")

    parser.add_argument(
        "-t",
        "--num-threads",
        dest="num_threads",
        type=int,
        default=1,
        help="Number of parallel jobs to run (passed to joblib.Parallel; -1 uses all cores)",
    )

    return parser.parse_args()

def run_nerpa_for_bgc_and_its_nrps(
        nerpa_root: Path,
        nerpa_script: Path,
        bgc_id: str,
        nrp_ids: List[str],
        pnrpdb: pl.DataFrame,
        antismash_results_all: Path,
        main_output_dir: Path,
        idx: int,
        total: int
) -> None:
        print(f'{idx}/{total}. Matching BGC {bgc_id} against {len(nrp_ids)} NRPs:\n{nrp_ids}')

        # Write the relevant NRPs to a temporary TSV file
        smiles_tsv_path: Path = nerpa_root / 'tmp' / f'{bgc_id}_nrps.tsv'
        smiles_tsv_path.parent.mkdir(parents=True, exist_ok=True)
        nrps_table: pl.DataFrame = (
            pnrpdb.filter(pl.col('ID').is_in(nrp_ids))
        )
        nrps_table.write_csv(smiles_tsv_path, separator='\t')

        antismash_results = antismash_results_all / bgc_id
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
            "--draw-hmms"
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

        subprocess.run(command, check=True)
    


def main():
    nerpa_root = Path(__file__).resolve().parent.parent
    assert nerpa_root.name.startswith('nerpa'), f"Expected nerpa_dir to be the nerpa repo, got {nerpa_root}"
    nerpa_script = nerpa_root / "nerpa.py"


    args = parse_args(nerpa_root)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    mibig_bgcs_info = MIBiG_BGCs_Info.load_from_nerpa_root(nerpa_root)
    pnrpdb_info = PNRPDB_Info.load_from_nerpa_root(nerpa_root)
    good_bgcs_to_nrps: Dict[str, List[str]] = get_good_mibig_bgs_to_nrps(
        mibig_bgcs_info,
        pnrpdb_info,
        min_num_a_domains=args.min_a_domains,
        min_num_nrps_monomers=args.min_monomers
    )

    pnrpdb = pl.read_csv(nerpa_root / 'data/input/pnrpdb2_expanded.tsv', separator='\t')
    local_paths = yaml.safe_load((nerpa_root / "local_paths.yaml").read_text())
    antismash_results_all = Path(local_paths['as_results_mibig4_nrps'])

    Parallel(n_jobs=args.num_threads)(
        delayed(run_nerpa_for_bgc_and_its_nrps)(
            nerpa_root=nerpa_root,
            nerpa_script=nerpa_script,
            bgc_id=bgc_id,
            nrp_ids=nrp_ids,
            pnrpdb=pnrpdb,
            antismash_results_all=antismash_results_all,
            main_output_dir=args.output_dir,
            idx=i,
            total=len(good_bgcs_to_nrps)
        )
        for i, (bgc_id, nrp_ids) in enumerate(good_bgcs_to_nrps.items(), start=1)
    )


if __name__ == '__main__':
    main()
