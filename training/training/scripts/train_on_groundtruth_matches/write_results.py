from typing import List
from pathlib import Path
import yaml
import pandas as pd
from prettytable import PrettyTable
from calculate_parameters import TrainedParameters
from dataclasses import asdict


def show_match(match: dict) -> str:
    result = ''
    result += f"Genome: {match['Genome']}\n"
    result += f"BGC_variant_idx: {match['BGC_variant_idx']}\n"
    result += f"NRP: {match['NRP']}\n"
    result += f"NRP_variant_idx: {match['NRP_variant_idx']}\n"
    result += f"NormalisedScore: {match['NormalisedScore']}\n"
    result += f"Score: {match['Score']}\n"
    result += f"Alignment:\n"
    for i, alignment in enumerate(match['Alignments']):
        if len(match['Alignments']) > 1:
            result += f'Fragment_{i}\n'
        fields = ['Gene', 'A-domain_idx', 'Top_scoring_residues', 'Modifying_domains', 'NRP_residue', 'NRP_chirality',
                  'NRP_modifications', 'rBAN_name', 'rBAN_idx', 'Alignment_step',
                  'Score', 'ResidueScore', 'MethylationScore', 'ChiralityScore']
        t = PrettyTable(fields, align='l', border=False)
        t.add_rows([alignment_step[field] for field in fields] for alignment_step in alignment)
        result += str(t) + '\n'

    return result


def write_results(matches: List[dict],
                  matches_table: pd.DataFrame,
                  parameters: TrainedParameters,
                  output_dir: Path):
    nrp_ids_good_matches = matches_table[matches_table['Verdict'].isin('good', 'was corrected')]['NRP']
    nrp_ids_to_invesigate = matches_table[matches_table['Verdict'] == 'for investigation']['NRP']
    nrp_ids_to_correct = matches_table[matches_table['Verdict'] == 'to be corrected']['NRP']

    good_matches = [match for match in matches if match['NRP'] in nrp_ids_good_matches]

    output_dir.mkdir(exist_ok=True, parents=True)
    # q: write good matches in a yaml file
    with open(output_dir / 'good_matches.yaml', 'w') as f:
        yaml.dump(good_matches, f)

    # q: write "to investigate" matches in a human-readable format
    with open(output_dir / 'matches_to_investigate.txt', 'w') as f:
        for match in matches:
            if match['NRP'] in nrp_ids_to_invesigate:
                f.write(show_match(match) + '\n')

    # q: write "to correct" matches in a human-readable format
    with open(output_dir / 'matches_to_correct.txt', 'w') as f:
        for match in matches:
            if match['NRP'] in nrp_ids_to_correct:
                f.write(show_match(match) + '\n')

    # q: write parameters in a yaml file
    with open(output_dir / 'parameters.yaml', 'w') as f:
        yaml.dump(asdict(parameters), f)