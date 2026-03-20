from pathlib import Path

import yaml

from src.antismash_parsing.bgc_variant_types import BGC_Variant
from src.rban_parsing.rban_parser import Parsed_rBAN_Record
from src.rban_parsing.retrieve_nrp_variants import rban_records_to_nrp_variants
from src.testing.simplified_alignment import simplified_alignment_to_light_alignments
from src.testing.testing_types import TestMatch
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Test simplified_alignment_to_light_alignments function.')
    parser.add_argument('nrp_id', type=str, help='The NRP ID to test, e.g., BGC0000374.6')
    args = parser.parse_args()
    if not args.nrp_id.startswith('BGC'):
        args.nrp_id = 'BGC' + '0' * (12 - len('BGC') - len(args.nrp_id) ) + args.nrp_id

    return args


def run_nerpa_single_mibig_nrp(nrp_id: str, nerpa_dir: Path):
    script = nerpa_dir / "_run_configurations/_nerpa_single_mibig_nrp.py"
    command = ["python3", str(script), nrp_id, '--dont-open-html']
    print(f'Running command: {" ".join(command)}')
    import subprocess
    subprocess.run(command, check=True)


def main():
    # TODO: accept nrp_id as a CL argument and run _nerpa_single_mibig_nrp.py with it at the start
    nerpa_dir = Path(__file__).parent.parent
    assert nerpa_dir.name.startswith(
        'nerpa'), f'Wrong nerpa_dir: {nerpa_dir}'
    args = parse_args()
    nrp_id = args.nrp_id

    approved_matches_yaml = nerpa_dir / 'data/for_training_and_testing/approved_matches.yaml'
    with approved_matches_yaml.open("r", encoding="utf-8") as f:
        approved_matches: list[TestMatch] = [TestMatch.from_yaml_dict(record)
                                             for record in yaml.safe_load(f)]

    test_match = next((match for match in approved_matches
                       if match.nrp_id == nrp_id), None)

    run_nerpa_single_mibig_nrp(nrp_id, nerpa_dir)
    nerpa_output_dir = nerpa_dir / 'nerpa_results/trash_single'
    bgc_variants_path = (nerpa_output_dir
                         / 'preprocessed_input'
                         / 'BGC_variants.yaml')
    with open(bgc_variants_path, 'r') as f:
        bgc_variants = [BGC_Variant.from_yaml_dict(record)
                        for record in yaml.safe_load(f)]

    bgc_id = nrp_id.split('.')[0]
    bgc_variant = next((variant for variant in bgc_variants
                        if variant.bgc_variant_id._genome_id() == bgc_id), None)

    # if bgc_variant is None:
    #     raise ValueError(f'BGC variant with ID {bgc_id} not found in {bgc_variants_path}')

    parsed_rban_records_yaml = (nerpa_output_dir
                         / 'preprocessed_input'
                         / 'parsed_rban_records.yaml')
    with open(parsed_rban_records_yaml, 'r') as f:
        parsed_rban_records = [Parsed_rBAN_Record.from_dict(record)
                               for record in yaml.safe_load(f)]

    nrp_variants = rban_records_to_nrp_variants(parsed_rban_records)
    nrp_variant = next((variant for variant in nrp_variants
                        if variant.nrp_variant_id.nrp_id == nrp_id), None)
    if nrp_variant is None:
        raise ValueError(f'NRP variant with ID {nrp_id} not found in {parsed_rban_records_yaml}')
    simplified_alignment = test_match.true_alignment
    for bgc_variant in bgc_variants:
        print(f'Trying BGC variant {bgc_variant.bgc_variant_id._genome_id()}-{bgc_variant.bgc_variant_id.variant_idx}...')
        light_alignments = list(simplified_alignment_to_light_alignments(
            simplified_alignment,
            bgc_variant,
            nrp_variant,
            debug=True
        ))
        if light_alignments:
            print(f'Found {len(light_alignments)} light alignments')
            # print('First one:')
            # for step in light_alignments[0]:
            #     bgc_module_str = f'Module {step.bgc_module_idx}' if step.bgc_module_idx is not None else 'No BGC module'
            #     nrp_monomer_str = f'Monomer {step.nrp_monomer_idx} ({step.nrp_monomer_residue})' if step.nrp_monomer_idx is not None else 'No NRP monomer'
            #     print(f'  {bgc_module_str} <-> {nrp_monomer_str}')
        else:
            print('No light alignments found')


if __name__ == "__main__":
    main()
