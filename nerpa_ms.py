import shutil
from typing import List, Literal
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
from src.monomer_names_helper import MonomerNamesHelper
from src.config import load_monomer_names_helper, load_config
import subprocess
import yaml

def parse_args() -> argparse.Namespace:
    # Detect if nerpa_results are provided before constructing the complete argument parser.
    prelim_parser = argparse.ArgumentParser(add_help=False)
    prelim_parser.add_argument("--nerpa-results", type=Path)
    prelim_args, _unknown_args = prelim_parser.parse_known_args()
    use_existing_results: bool = prelim_args.nerpa_results is not None

    parent_parsers: list[argparse.ArgumentParser] = []
    if not use_existing_results:
        nerpa_default_config = load_nerpa_config()
        parent_parsers.append(
            build_nerpa_args_parser(nerpa_default_config, add_help=False)
        )

    nerpa_ms_parser = argparse.ArgumentParser(
        description="Run NERPA-MS pipeline",
        parents=parent_parsers,
        add_help=True,
    )

    if use_existing_results:
        nerpa_ms_parser.add_argument(
            "--output-dir",
            type=Path,
            required=True,
            help="Directory in which to write NERPA-MS results",
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
        help="Maximum number of Nerpa matches to consider for candidate generation (default: 10)"
    )
    nerpa_ms_parser.add_argument(
        '--nerpa-results',
        type=Path,
        help="Path to the directory with Nerpa results"
    )
    nerpa_ms_parser.add_argument(
        '--nerpa-match-ids-for-candidate-generation',
        type=str,
        help=(
            'The 0-based match IDs used for candidate generation in the format of a list of indices and ranges (e.g. "0, 2-4, 5"). '
            'Overrides --max-nerpa-matches-for-candidate-generation. '
            'Requires --nerpa-results'
        )
    )
    nerpa_ms_parser.add_argument(
        '--max-spectra-matches-per-nerpa-match',
        type=int,
        help="Maximum number of spectra matches to display per Nerpa match in the HTML report (default: 5)",
        default=5
    )
    nerpa_ms_parser.add_argument(
        '--kakapo-mode',
        help=(
            "Kakapo modification matching mode (default: free-mod). "
            "Choices: no-mod (unmodified only), "
            "free-mod (single modification on any monomer), "
            "x-mod (single modification only on X monomers), "
            "two-mod (known modification on base monomer + residual on X monomer)."
        ),
        choices=['no-mod', 'free-mod', 'x-mod', 'two-mod'],
        default='free-mod'
    )

    args = nerpa_ms_parser.parse_args()

    if args.nerpa_results is None and args.nerpa_match_ids_for_candidate_generation is not None:
        nerpa_ms_parser.error(
            "--nerpa-match-ids-for-candidate-generation requires --nerpa-results"
        )

    args.output_dir = args.output_dir.resolve()
    args.spectra = args.spectra.resolve()

    return args


def nerpa_ms_core(
        nerpa_root: Path,
        nerpa_results: Path,
        output_dir: Path,
        nerpa_match_ids: str | None = None,
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
    if nerpa_match_ids is not None:
        command.extend(['--nerpa-match-ids', nerpa_match_ids])

    print(f"Running in {nerpa_ms_root}:\n  {' '.join(command)}")
    subprocess.run(command, cwd=nerpa_ms_root, check=True)


def nerpa_ms_varquest(
        nerpa_root: Path,
        nerpa_ms_core_results_dir: Path,
        spectra_path: Path,
        output_dir: Path,
        mode: Literal['no-mod', 'free-mod', 'x-mod', 'two-mod'] = 'free-mod'
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
        '-f', 'json',
        '-m', mode
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
        nerpa_ms_varquest_results: Path,
        monomer_names_helper: MonomerNamesHelper,
        max_spectra_matches_per_nerpa_match: int,
) -> None:
   cfg = HTMLReportConfig(
        main_out_dir=nerpa_results,
        nerpa_root=nerpa_root,
        mode='nerpa-ms',
        generated_candidate_nrps_path=nerpa_ms_core_results / "nerpa-ms-output.json",
        mass_spec_matching_results_dir=nerpa_ms_varquest_results,
        max_spectra_matches_per_nerpa_match=max_spectra_matches_per_nerpa_match,
    )
   bgc_variants_info = load_bgc_variants_info(nerpa_results)
   nrp_variants_info = load_nrp_variants_info(nerpa_results)
   matches = load_matches(nerpa_results)

   create_html_report(
        bgc_variants_info=bgc_variants_info,
        nrp_variants_info=nrp_variants_info,
        matches=matches,
        cfg=cfg,
        monomer_names_helper=monomer_names_helper,
        debug_output=False,
    )
    

def get_monomer_names_helper() -> MonomerNamesHelper:
    cfg = load_config()
    return load_monomer_names_helper(
        cfg.monomers_config,
        cfg.nerpa_dir
    )
       
def main():
    nerpa_root = Path(__file__).parent.resolve()
    args = parse_args()

    if args.nerpa_results is None:
        nerpa_main(args=args)
    else:
        if args.output_dir.exists():
            print(f"Removing existing output directory: {args.output_dir}")
            shutil.rmtree(args.output_dir)

        print(f"Copying Nerpa results from: {args.nerpa_results}")
        shutil.copytree(args.nerpa_results, args.output_dir)

    nerpa_ms_core_out = args.output_dir / "nerpa_ms_results"
    nerpa_ms_core(
        nerpa_root=nerpa_root,
        nerpa_results=args.output_dir,
        output_dir=nerpa_ms_core_out,
        max_edits=args.max_edits,
        num_variants_per_num_edits=args.num_variants_per_num_edits,
        max_nerpa_matches=args.max_nerpa_matches_for_candidate_generation,
        nerpa_match_ids=args.nerpa_match_ids_for_candidate_generation,
    )

    nerpa_ms_varquest_out = args.output_dir / "nerpa_ms_varquest_results"
    nerpa_ms_varquest(
        nerpa_root=nerpa_root,
        nerpa_ms_core_results_dir=nerpa_ms_core_out,
        spectra_path=args.spectra,
        output_dir=nerpa_ms_varquest_out,
        mode=args.kakapo_mode
    )
    print(f"Nerpa-MS results written to: {nerpa_ms_varquest_out}")

    monomer_names_helper = get_monomer_names_helper()
    build_html_report(
        nerpa_root=nerpa_root,
        nerpa_results=args.output_dir,
        nerpa_ms_core_results=nerpa_ms_core_out,
        nerpa_ms_varquest_results=nerpa_ms_varquest_out,
        monomer_names_helper=monomer_names_helper,
        max_spectra_matches_per_nerpa_match=args.max_spectra_matches_per_nerpa_match
    )

if __name__ == "__main__":
    main()
