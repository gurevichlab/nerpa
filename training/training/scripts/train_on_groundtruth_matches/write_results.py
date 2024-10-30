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
from src.matching.matching_types_match import Match
from src.data_types import BGC_Variant, NRP_Variant


def get_fragmented_matches_nrp_ids(matches: List[Match],
                                   bgc_variants_by_nrp_id: Dict[str, BGC_Variant],
                                   nrp_variants_by_nrp_id: Dict[str, NRP_Variant]) -> List[str]:
    nrp_ids = {match.nrp_variant_info.nrp_id for match in matches}
    def bgc_fragments_cnt(nrp_id: str) -> int:
        return max(module.fragment_idx
                   for module in bgc_variants_by_nrp_id[nrp_id].modules) + 1
    def nrp_fragments_cnt(nrp_id: str) -> int:
        return len(nrp_variants_by_nrp_id[nrp_id].fragments)
    return [nrp_id for nrp_id in nrp_ids
            if bgc_fragments_cnt(nrp_id) > 1
            or nrp_fragments_cnt(nrp_id) > 1]


def write_results(matches: List[Match],
                  bgc_variants_by_nrp_id: Dict[str, BGC_Variant],
                  nrp_variants_by_nrp_id: Dict[str, NRP_Variant],
                  matches_table: pd.DataFrame,
                  parameters: TrainedParameters,
                  output_dir: Path):
    nrp_ids_good_matches = list(matches_table[matches_table['Verdict'].isin(['good', 'was corrected'])]['NRP variant'])
    nrp_ids_to_investigate = list(matches_table[matches_table['Verdict'] == 'for investigation']['NRP variant'])
    nrp_ids_to_correct = list(matches_table[matches_table['Verdict'] == 'to be corrected']['NRP variant'])
    nrp_ids_fragmented_matches = get_fragmented_matches_nrp_ids(matches,
                                                                bgc_variants_by_nrp_id,
                                                                nrp_variants_by_nrp_id)

    nrp_id_to_match = {match.nrp_variant_info.nrp_id: match for match in matches}

    nrp_ids_good_matches = set(nrp_ids_good_matches) & set(nrp_id_to_match)
    nrp_ids_to_investigate = set(nrp_ids_to_investigate) & set(nrp_id_to_match)
    nrp_ids_to_correct = set(nrp_ids_to_correct) & set(nrp_id_to_match)
    nrp_ids_fragmented_matches = set(nrp_ids_fragmented_matches) & set(nrp_id_to_match)

    nrp_ids_to_investigate_many_fragments = (set(nrp_ids_to_investigate) | set(nrp_ids_to_correct)) & set(nrp_ids_fragmented_matches)
    good_matches = [nrp_id_to_match[nrp_id] for nrp_id in nrp_ids_good_matches
                    if nrp_id in nrp_id_to_match]

    output_dir.mkdir(exist_ok=True, parents=True)
    # q: write good matches in a yaml file
    write_yaml(good_matches, output_dir / 'good_matches.yaml')

    # q: write "to investigate" matches in a human-readable format
    with open(output_dir / 'matches_to_investigate.txt', 'w') as f:
            for nrp_id in nrp_ids_to_investigate:
                f.write(str(nrp_id_to_match[nrp_id]) + '\n')

    # q: write "to correct" matches in a human-readable format
    with open(output_dir / 'matches_to_correct.txt', 'w') as f:
            for nrp_id in nrp_ids_to_correct:
                f.write(str(nrp_id_to_match[nrp_id]) + '\n')

    with open(output_dir / 'matches_to_investigate_many_fragments.txt', 'w') as f:
        for nrp_id in nrp_ids_to_investigate_many_fragments:
            f.write(str(nrp_id_to_match[nrp_id]) + '\n')

    nrp_ids_list = list(nrp_ids_to_investigate_many_fragments)
    block_size = (len(nrp_ids_to_investigate_many_fragments) + 3) // 4  # round up
    for i, name in enumerate(['Azat', 'Ilia', 'Lesha', 'Sasha']):
        nrp_ids_for_person = nrp_ids_list[i * block_size: (i + 1) * block_size]
        with open(output_dir / f'matches_to_investigate_{name}_nrp_ids.txt', 'w') as f:
            f.write('\n'.join(nrp_ids_for_person))
        with open(output_dir / f'matches_to_investigate_{name}_alignments.txt', 'w') as f:
            for nrp_id in nrp_ids_for_person:
                f.write(str(nrp_id_to_match[nrp_id]) + '\n')


    # q: write parameters in a yaml file
    write_yaml(asdict(parameters), output_dir / 'parameters.yaml')  # asdict does not work for nested dataclasses
