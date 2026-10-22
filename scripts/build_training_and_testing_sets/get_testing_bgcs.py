from collections import defaultdict

from pandas.compat.pyarrow import pa
from src.benchmarking.data_frames import (
    MIBiG_BGCs_Info,
    PNRPDB_Info,
)
import polars as pl
from pathlib import Path
from src.testing.testing_types import CuratedAlignment
import yaml
import argparse

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Write antiSMASH paths for BGCs used for testing:\n"
            "MIBiG BGCs dissimilar to the BGCs used in training, with enough A domains and their corresponding NRPs with enough NRPS monomers."
        )
    )
    parser.add_argument("--curated-alignments", type=Path)
    parser.add_argument("--antismash-results", type=Path)
    parser.add_argument("--pnrpdb-info", type=Path)
    parser.add_argument("--mibig-bgcs-info", type=Path)
    parser.add_argument("--min-a-domains", type=int, default=3)
    parser.add_argument("--min-nrps-monomers", type=int, default=3)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main():
    args = parse_arguments()
    
    mibig_bgcs_info = MIBiG_BGCs_Info.from_csv(args.mibig_bgcs_info)
    curated_alignments: list[CuratedAlignment] = [
        CuratedAlignment.from_yaml_dict(test_match_dict)
        for test_match_dict in yaml.safe_load(args.curated_alignments.open('r'))
    ]

    training_bgcs: list[str] = sorted(ca.mibig_id for ca in curated_alignments)

    training_bgcs_info: pl.DataFrame = mibig_bgcs_info.filter(
        pl.col(MIBiG_BGCs_Info.BGC_ID).is_in(training_bgcs)
    )
    training_bgcs_fams: set[str] = set(training_bgcs_info[MIBiG_BGCs_Info.BIGSCAPE_FAMILIES])

    dissimilar_to_training_bgcs_info: pl.DataFrame = (
        mibig_bgcs_info
        .filter(~pl.col(MIBiG_BGCs_Info.BIGSCAPE_FAMILIES).is_in(training_bgcs_fams))
    )

    pnrpdb_info: PNRPDB_Info = PNRPDB_Info.from_csv(args.pnrpdb_info)
    # MIBiG or NORINE NRPs with at least 3 NRPS monomers
    nrps_with_enough_mons: set[str] = set(
        pnrpdb_info
        .filter(pl.col(PNRPDB_Info.NUM_NRPS_MONOMERS) >= args.min_nrps_monomers)
        [PNRPDB_Info.COMPOUND_ID]
    )

    test_bgcs: set[str] = set(
        dissimilar_to_training_bgcs_info
        .filter(pl.col(MIBiG_BGCs_Info.NUM_A_DOMAINS) >= args.min_a_domains)
        [MIBiG_BGCs_Info.BGC_ID]
    )

    bgc_to_nrps: dict[str, set[str]] = defaultdict(set)
    for row in mibig_bgcs_info.iter_rows(named=True):
        bgc_id = row[MIBiG_BGCs_Info.BGC_ID]
        nrp_ids = (
            row[MIBiG_BGCs_Info.NRP_IDS].split(',')
            if row[MIBiG_BGCs_Info.NRP_IDS]
            else []
        )
        bgc_to_nrps[bgc_id] = set(nrp_ids)

    # Keep only BGCs that have at least one NRP with enough NRPS monomers
    test_bgcs: set[str] = {
        bgc
        for bgc in test_bgcs
        if any(nrp in nrps_with_enough_mons
               for nrp in bgc_to_nrps[bgc])
    }

    test_bgcs -= {'BGC0001561'}  # a mistake in MIBiG (TODO: elaborate)

    individual_paths: list[Path] = [
        args.antismash_results / bgc_id
        for bgc_id in sorted(test_bgcs)
    ]

    output_text = "\n".join(map(str, individual_paths))
    if output_text:
        output_text += "\n"

    args.output.write_text(output_text, encoding="utf-8")

    readme_path = Path(f"{args.output}.readme")
    readme_text = (
        "This file contains one antiSMASH result path per line.\n\n"
        "Included BGCs:\n"
        "- are not in the same BiG-SCAPE families as the training BGCs;\n"
        f"- contain at least {args.min_a_domains} A domains; and\n"
        "- have at least one associated NRP containing "
        f"at least {args.min_nrps_monomers} NRPS monomers.\n\n"
        f"Number of BGCs: {len(test_bgcs)}\n"
    )
    readme_path.write_text(readme_text, encoding="utf-8")

    
if __name__ == "__main__":
    main()
