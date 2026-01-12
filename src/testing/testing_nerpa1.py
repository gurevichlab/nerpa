from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List

from src.antismash_parsing.antismash_parser_types import GeneId
from src.antismash_parsing.bgc_variant_types import A_Domain_ID
from src.monomer_names_helper import NorineMonomerName
from src.testing.simplified_alignment import (
    SimplifiedAlignment,
)
from src.testing.simplified_match import SimplifiedMatch


def nerpa1_match_to_simplified_match(nerpa1_match: str) -> SimplifiedMatch:
    lines = nerpa1_match.strip().split('\n')
    nrp_id = lines[0].strip().split('_')[0]  # ignore variant idx
    bgc_path = lines[1].strip()
    bgc_id = Path(bgc_path).parts[-1].split('_')[0]  # Extract BGC ID from path

    # Find the start of the ALIGNMENT data (skip metadata and header)
    alignment_start = next(i for i, line in enumerate(lines) if line.strip() == 'ALIGNMENT:') + 1  # +1 to get to header
    header_line = lines[alignment_start]
    header = header_line.split()
    expected_header = [
        'ORF_ID',
        'A_domain_Idx',
        'Prediction_DL-config',
        'Prediction_Top_Residue',
        'Prediction_Top_Score',
        'Prediction_Modifications',
        'Matched_Residue',
        'Matched_Residue_Score',
        'Nerpa_Score',
        'Monomer_Idx',
        'Monomer_Code',
        'Monomer_DL-config',
        'Monomer_Residue',
        'Monomer_Modifications',
   ]
    if header != expected_header:
        raise ValueError(f"Unexpected header: {header}")

    data_lines = lines[alignment_start + 1:]

    simplified_alignment = []
    for line in data_lines:
        if not line.strip():
            continue
        step = dict(zip(header, line.split()))
        if step['ORF_ID'] != '-':
            gene_id = step['ORF_ID'].split('_', 1)[1]  # remove 'ctgX_' prefix
            a_domain_idx = int(step['A_domain_Idx'])
            a_domain_id = A_Domain_ID(GeneId(gene_id), a_domain_idx)
        else:
            a_domain_id = None

        if (mon_code := step['Monomer_Code']) != '-':
            if mon_code.startswith('@L-') or mon_code.startswith('@D-'):  # restore original rBAN names
                mon_code = mon_code[3:]
            monomer_name = NorineMonomerName(mon_code)
        else:
            monomer_name = None

        simplified_alignment.append((a_domain_id, monomer_name))

    return SimplifiedMatch(bgc_id=bgc_id,
                           nrp_id=nrp_id,
                           alignment=simplified_alignment)


def load_nerpa1_matches(nerpa_results_dir: Path) -> List[SimplifiedMatch]:
    matches = []
    for match_file in (nerpa_results_dir / 'details').glob('*.match'):
        content = match_file.read_text()
        match_strs = [mstr for mstr in content.split('\n\n') if mstr.strip()]
        matches.extend(nerpa1_match_to_simplified_match(match_str)
                       for match_str in match_strs)
    
    return matches


def test_from_nerpa1_match():
    # Path to the provided data file
    data_file = Path('/home/ilianolhin/git/nerpa2/data/for_training_and_testing/nerpa1_results/details/BGC0000290.0_variant0.match')
    with open(data_file, 'r') as f:
        content = f.read()
    print(f'Parsing the first match from {data_file}')

    # Split into blocks (assuming blocks are separated by double newlines)
    blocks = content.strip().split('\n\n')
    first_block = blocks[0]
    
    # Parse the first block
    match = nerpa1_match_to_simplified_match(first_block)
    print(f'Parsed match: {match}')

    # Expected alignment for the first block
    expected_alignment = [
        (A_Domain_ID(GeneId('teiA'), 0), None),
        (A_Domain_ID(GeneId('teiA'), 1), NorineMonomerName('Cl-Tyr')),
        (A_Domain_ID(GeneId('teiB'), 0), NorineMonomerName('Dhpg')),
        (A_Domain_ID(GeneId('teiC'), 0), NorineMonomerName('Hpg')),
        (A_Domain_ID(GeneId('teiC'), 1), NorineMonomerName('Cl-Hpg')),
        (A_Domain_ID(GeneId('teiC'), 2), NorineMonomerName('bOH-Cl-Tyr')),
        (A_Domain_ID(GeneId('teiD'), 0), NorineMonomerName('Dhpg')),
    ]
    
    if all([match.nrp_id == 'BGC0000290.0',
            match.bgc_id == 'BGC0000441',
            match.alignment == expected_alignment]):
        print('Match parsed correctly!')
    else:
        print('Mismatch in parsed match!')
        print(f'Expected nrp_id: BGC0000290.0, got: {match.nrp_id}')
        print(f'Expected bgc_id: BGC0000441, got: {match.bgc_id}')
        print(f'Expected alignment: {expected_alignment}, got: {match.alignment}')



if __name__ == '__main__':
    test_from_nerpa1_match()