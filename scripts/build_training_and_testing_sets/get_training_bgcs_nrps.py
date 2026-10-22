"""Write antiSMASH result paths for BGCs listed in approved matches."""

import argparse
from pathlib import Path
import polars as pl
from src.testing.testing_types import CuratedAlignment

import json


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write antiSMASH paths for approved BGC matches."
    )
    parser.add_argument("--curated-alignments", type=Path)
    parser.add_argument("--pnrpdb-expanded", type=Path)
    parser.add_argument("--antismash-results", type=Path)
    parser.add_argument("--output-bgcs", type=Path)
    parser.add_argument("--output-nrps", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()

    curated_alignments_data = json.loads(args.curated_alignments.read_text(encoding="utf-8"))
    curated_alignments: list[CuratedAlignment] = [
        CuratedAlignment.from_yaml_dict(test_match_dict)
        for test_match_dict in curated_alignments_data["curated_alignments"]
    ]

    training_bgcs: list[str] = sorted({ca.mibig_id for ca in curated_alignments})
    individual_paths: list[Path] = [
        args.antismash_results / bgc_id
        for bgc_id in training_bgcs
    ]

    training_nrps: set[str] = {ca.nrp_id for ca in curated_alignments}

    bgcs_paths_text = "\n".join(map(str, individual_paths))
    args.output_bgcs.write_text(bgcs_paths_text, encoding="utf-8")

    bgcs_readme_path = Path(f"{args.output_bgcs}.readme")
    readme_text = (
        "This file lists antiSMASH result directories for BGCs found in the "
        "curated alignments.\n"
        "Each line contains one directory path.\n"
        f"Curated alignments: {args.curated_alignments}\n"
        f"antiSMASH results: {args.antismash_results}\n"
    )
    bgcs_readme_path.write_text(readme_text, encoding="utf-8")

    pnrpdb_expanded: pl.DataFrame = pl.read_csv(args.pnrpdb_expanded, separator='\t')
    pnrpdb_filtered: pl.DataFrame = (
        pnrpdb_expanded
        .filter(pl.col("ID").is_in(training_nrps))
    )

    pnrpdb_filtered.write_csv(args.output_nrps, separator='\t')
    nrps_readme_path = Path(f"{args.output_nrps}.readme")
    nrps_readme_text = (
        "This file lists NRPs found in the curated alignments.\n"
        "Each line contains one NRP ID.\n"
        f"Curated alignments: {args.curated_alignments}\n"
        f"PNRPDB expanded: {args.pnrpdb_expanded}\n"
    )
    nrps_readme_path.write_text(nrps_readme_text, encoding="utf-8")
    
        
if __name__ == "__main__":
    main()

