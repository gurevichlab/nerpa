from typing import List
import argparse
from pathlib import Path
from nerpa import main as nerpa_main
from src.config import load_config as load_nerpa_config
from src.pipeline.command_line_args_helper import build_cmdline_args_parser as build_nerpa_args_parser
from src.rban_parsing.rban_parser import Parsed_rBAN_Record
from src.antismash_parsing.bgc_variant_types import (
    BGC_Variant,
    BGC_Variants_Info,
    BGC_Variant_ID,
)
from src.rban_parsing.nrp_variant_types import (
    NRP_Variant,
    NRP_Variants_Info,
    NRP_Variant_ID,
)
from src.build_output.html_reporter import HTMLReportConfig, create_html_report
from src.matching.match_type import Match
import subprocess
import yaml

def parse_args() -> argparse.Namespace:
    nerpa_default_config = load_nerpa_config()
    nerpa_parser: argparse.ArgumentParser = \
        build_nerpa_args_parser(nerpa_default_config, add_help=False)

    nerpa_ms_parser = argparse.ArgumentParser(
        description="Run NERPA-MS pipeline",
        parents=[nerpa_parser],
        add_help=True
    )
    nerpa_ms_parser.add_argument(
        '--spectra',
        type=Path,
        required=True,
        help="Path to the spectra file in MGF format"
    )
    nerpa_ms_parser.add_argument(
        '--max-edits',
        type=int,
        default=2,
        help="Maximum number of edits allowed in generating candidate NRP variants (default: 2)"
    )
    nerpa_ms_parser.add_argument(
        '--num-variants-per-num-edits',
        type=int,
        default=100,
        help="Maximum number of candidate NRP variants to generate per number of edits (default: 100)"
    )
    nerpa_ms_parser.add_argument(
        '--max-nerpa-matches-for-candidate-generation',
        type=int,
        default=10,
        help="Maximum number of NERPA matches to consider for candidate generation (default: 10)"
    )

    args = nerpa_ms_parser.parse_args()

    args.output_dir = args.output_dir.resolve()
    args.spectra = args.spectra.resolve()

    return args


def nerpa_ms_core(
        nerpa_root: Path,
        nerpa_results: Path,
        output_dir: Path,
        max_nerpa_matches: int = 100,
        max_edits: int = 2,
        num_variants_per_num_edits: int = 5
) -> None:
    nerpa_ms_root = (
        nerpa_root
        / 'src'
        / 'nerpa_ms'
        / 'nerpa_ms_core'
    )

    subprocess.run(["cargo", "build"], cwd=nerpa_ms_root, check=True)
    command = [
        'cargo',  'run',
        '--bin', 'nerpa_ms_core', '--',
        '--nerpa-results', str(nerpa_results),
        '--max-nerpa-matches', str(max_nerpa_matches),
        '--max-edits', str(max_edits),
        '--num-variants-per-num-edits', str(num_variants_per_num_edits),
        '--out', str(output_dir),
        '--monomers-db-json', str(nerpa_ms_root / 'data/monomers_db.json'),
        '--nerpa-root', str(nerpa_root),
    ]
    print(f"Running in {nerpa_ms_root}:\n  {' '.join(command)}")
    subprocess.run(command, cwd=nerpa_ms_root, check=True)


def nerpa_ms_varquest(
        nerpa_root: Path,
        nerpa_ms_core_results_dir: Path,
        spectra_path: Path,
        output_dir: Path
) -> None:
    # with open(nerpa_ms_core_results_dir / "nerpa-ms-output.json", 'r') as f:
    #     nerpa_ms_core_results = json.load(f)

    # generated_candidates: List[dict] = [
    #     candidate
    #     for _nerpa_match_id, candidates_for_match in nerpa_ms_core_results.items()
    #     for candidate in candidates_for_match['new_variants'].values()
    # ]

    # generated_candidates_path = nerpa_ms_core_results_dir / "nerpa-ms-output--flattened.json"
    # with open(generated_candidates_path, 'w') as f:
    #     json.dump(
    #         generated_candidates,
    #         f,
    #         indent=2
    #     )

    kakapo_exec: Path = nerpa_root / "envs" / "kakapo" / "bin" / "kakapo"
    generated_candidates_path = nerpa_ms_core_results_dir / "nerpa-ms-output.json"

    command = [
        str(kakapo_exec),
        '-s', str(spectra_path),
        '-str', str(generated_candidates_path),
        '-o', str(output_dir),
        '-f', 'json'
    ]
    print(f"Running Kakapo:\n  {' '.join(command)}")
    subprocess.run(command, check=True)


def load_bgc_variants_info(nerpa_results: Path) -> BGC_Variants_Info:
    bgc_variants_path = (
        nerpa_results
        / "preprocessed_input"
        / "BGC_variants.yaml"
    )
    with open(bgc_variants_path, 'r') as f:
        bgc_variants = [BGC_Variant.from_dict(bgc_dict) for bgc_dict in yaml.safe_load(f)]

    bgc_representatives_path = (
        nerpa_results
        / "preprocessed_input"
        / "BGC_representatives.yaml"
    )
    with open(bgc_representatives_path, 'r') as f:
        bgc_representatives_data = yaml.safe_load(f)
    bgc_id_to_repr_id = {
        BGC_Variant_ID.from_dict(item['bgc_id']):
        BGC_Variant_ID.from_dict(item['representative_id'])
        for item in bgc_representatives_data
    }

    return BGC_Variants_Info(
        bgc_variants=bgc_variants,
        bgc_id_to_repr_id=bgc_id_to_repr_id
    )


def load_nrp_variants_info(nerpa_results: Path) -> NRP_Variants_Info:
    parsed_rban_records_path = (
        nerpa_results
        / "preprocessed_input"
        / "parsed_rban_records.yaml"
    )
    with open(parsed_rban_records_path, 'r') as f:
        parsed_rban_records = [Parsed_rBAN_Record.from_dict(record) for record in yaml.safe_load(f)]

    nrp_variants_path = (
        nerpa_results
        / "preprocessed_input"
        / "NRP_variants.yaml"
    )
    with open(nrp_variants_path, 'r') as f:
        nrp_variants = [NRP_Variant.from_dict(record) for record in yaml.safe_load(f)]

    nrp_representatives_path = (
        nerpa_results
        / "preprocessed_input"
        / "NRP_representatives.yaml"
    )
    with open(nrp_representatives_path, 'r') as f:
        nrp_representatives_data = yaml.safe_load(f)
    nrp_id_to_repr_id = {
        NRP_Variant_ID(*item['nrp_id']):
        NRP_Variant_ID(*item['representative_id'])
        for item in nrp_representatives_data
    }

    return NRP_Variants_Info(
        nrp_variants=nrp_variants,
        rban_records=parsed_rban_records,
        nrp_id_to_repr_id=nrp_id_to_repr_id
    )


def load_matches(nerpa_results: Path) -> list[Match]:
    matches_path = (
        nerpa_results
        / "matches_details"
        / "matches.yaml"
    )
    with open(matches_path, 'r') as f:
        matches_data = yaml.safe_load(f)
    matches = [Match.from_dict(match_dict) for match_dict in matches_data]
    return matches

        
def build_html_report(
        nerpa_root: Path,
        nerpa_results: Path,
        nerpa_ms_core_results: Path,
        nerpa_ms_varquest_results: Path
) -> None:
   cfg = HTMLReportConfig(
        main_out_dir=nerpa_results,
        nerpa_root=nerpa_root,
        mode='nerpa-ms',
        generated_candidate_nrps_path=nerpa_ms_core_results / "nerpa-ms-output.json",
        mass_spec_matching_results_dir=nerpa_ms_varquest_results
    )
   bgc_variants_info = load_bgc_variants_info(nerpa_results)
   nrp_variants_info = load_nrp_variants_info(nerpa_results)
   matches = load_matches(nerpa_results)

   create_html_report(
        bgc_variants_info=bgc_variants_info,
        nrp_variants_info=nrp_variants_info,
        matches=matches,
        cfg=cfg,
        debug_output=False
    )
    
       
def main():
    nerpa_root = Path(__file__).parent.resolve()
    args = parse_args()

    nerpa_main(args=args)

    nerpa_ms_core_out = args.output_dir / "nerpa_ms_results"
    nerpa_ms_core(
        nerpa_root=nerpa_root,
        nerpa_results=args.output_dir,
        output_dir=nerpa_ms_core_out,
        max_edits=args.max_edits,
        num_variants_per_num_edits=args.num_variants_per_num_edits,
        max_nerpa_matches=args.max_nerpa_matches_for_candidate_generation
    )

    nerpa_ms_varquest_out = args.output_dir / "nerpa_ms_varquest_results"
    nerpa_ms_varquest(
        nerpa_root=nerpa_root,
        nerpa_ms_core_results_dir=nerpa_ms_core_out,
        spectra_path=args.spectra,
        output_dir=nerpa_ms_varquest_out
    )

    build_html_report(
        nerpa_root=nerpa_root,
        nerpa_results=args.output_dir,
        nerpa_ms_core_results=nerpa_ms_core_out,
        nerpa_ms_varquest_results=nerpa_ms_varquest_out
    )

if __name__ == "__main__":
    main()
