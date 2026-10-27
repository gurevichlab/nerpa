from src.benchmarking.data_frames import (
    PNRPDB_Info,
)
import polars as pl
from pathlib import Path
import argparse

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Write TSV table with the NRPs used for testing:\n"
            "MIBiG or Norine NRPs with enough NRPS monomers."
        )
    )
    parser.add_argument("--pnrpdb-expanded", type=Path)
    parser.add_argument("--pnrpdb-info", type=Path)
    parser.add_argument("--min-nrps-monomers", type=int, default=3)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def is_bgc_or_norine_nrp() -> pl.Expr:
    return (pl.col(PNRPDB_Info.COMPOUND_ID).str.starts_with('BGC') |
            pl.col(PNRPDB_Info.COMPOUND_ID).str.starts_with('NOR'))


def main():
    args = parse_arguments()
    
    pnrpdb_info: PNRPDB_Info = PNRPDB_Info.from_csv(args.pnrpdb_info)
    # MIBiG or NORINE NRPs with at least 3 NRPS monomers
    test_nrps: set[str] = set(
        pnrpdb_info
        .filter(is_bgc_or_norine_nrp() &
                (pl.col(PNRPDB_Info.NUM_NRPS_MONOMERS) >= args.min_nrps_monomers))
        [PNRPDB_Info.COMPOUND_ID]
    )
    pnrpdb_expanded: pl.DataFrame = pl.read_csv(args.pnrpdb_expanded, separator='\t')
    test_nrps_table: pl.DataFrame = (
        pnrpdb_expanded
        .filter(pl.col("ID").is_in(test_nrps))
    )
    
    test_nrps_table.write_csv(args.output, separator='\t')

    readme_path: Path = Path(f"{args.output}.readme")
    readme_path.write_text(
        "TSV table containing MIBiG and NORINE NRPs with at least "
        f"{args.min_nrps_monomers} NRPS monomers.\n",
        encoding="utf-8",
    )

if __name__ == "__main__":
    main()
