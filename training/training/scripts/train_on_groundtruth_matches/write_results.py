from typing import Dict, List
from pathlib import Path
import yaml
import pandas as pd
from prettytable import PrettyTable
from calculate_parameters import TrainedParameters
from dataclasses import asdict
from src.write_results import write_yaml
from step_function import plot_step_function
from itertools import chain


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


def get_fragmented_matches_nrp_ids(matches: List[dict],
                                   bgc_variants: Dict[str, dict],
                                   nrp_variants: Dict[str, dict]) -> List[str]:
    nrp_ids = {match['NRP'] for match in matches}
    return [nrp_id for nrp_id in nrp_ids
            if len(bgc_variants[nrp_id]['fragments']) > 1
            or len(nrp_variants[nrp_id]['fragments']) > 1]


def write_results(matches: List[dict],
                  bgc_variants: Dict[str, dict],
                  nrp_variants: Dict[str, dict],
                  matches_table: pd.DataFrame,
                  parameters: TrainedParameters,
                  output_dir: Path):
    nrp_ids_good_matches = list(matches_table[matches_table['Verdict'].isin(['good', 'was corrected'])]['NRP variant'])
    nrp_ids_to_investigate = list(matches_table[matches_table['Verdict'] == 'for investigation']['NRP variant'])
    nrp_ids_to_correct = list(matches_table[matches_table['Verdict'] == 'to be corrected']['NRP variant'])
    nrp_ids_fragmented_matches = get_fragmented_matches_nrp_ids(matches,
                                                                bgc_variants,
                                                                nrp_variants)

    nrp_id_to_match = {match['NRP']: match for match in matches}
    nrp_ids_to_investigate_many_fragments = (set(nrp_ids_to_investigate) | set(nrp_ids_to_correct)) & set(nrp_ids_fragmented_matches)
    good_matches = [nrp_id_to_match[nrp_id] for nrp_id in nrp_ids_good_matches]

    output_dir.mkdir(exist_ok=True, parents=True)
    # q: write good matches in a yaml file
    write_yaml(good_matches, output_dir / 'good_matches.yaml')

    # q: write "to investigate" matches in a human-readable format
    with open(output_dir / 'matches_to_investigate.txt', 'w') as f:
            for nrp_id in nrp_ids_to_investigate:
                f.write(show_match(nrp_id_to_match[nrp_id]) + '\n')

    # q: write "to correct" matches in a human-readable format
    with open(output_dir / 'matches_to_correct.txt', 'w') as f:
            for nrp_id in nrp_ids_to_correct:
                f.write(show_match(nrp_id_to_match[nrp_id]) + '\n')

    with open(output_dir / 'matches_to_investigate_many_fragments.txt', 'w') as f:
        for nrp_id in nrp_ids_to_investigate_many_fragments:
            f.write(show_match(nrp_id_to_match[nrp_id]) + '\n')

    nrp_ids_list = list(nrp_ids_to_investigate_many_fragments)
    block_size = (len(nrp_ids_to_investigate_many_fragments) + 3) // 4  # round up
    for i, name in enumerate(['Azat', 'Ilia', 'Lesha', 'Sasha']):
        nrp_ids_for_person = nrp_ids_list[i * block_size: (i + 1) * block_size]
        with open(output_dir / f'matches_to_investigate_{name}_nrp_ids.txt', 'w') as f:
            f.write('\n'.join(nrp_ids_for_person))
        with open(output_dir / f'matches_to_investigate_{name}_alignments.txt', 'w') as f:
            for nrp_id in nrp_ids_for_person:
                f.write(show_match(nrp_id_to_match[nrp_id]) + '\n')


    # q: write parameters in a yaml file
    write_yaml(asdict(parameters), output_dir / 'parameters.yaml')
