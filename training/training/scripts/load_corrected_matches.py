# load approved matches
from csv import DictReader
from yaml import safe_load, safe_dump
from pathlib import Path

approved_matches_table = Path('/home/ilianolhin/git/nerpa2/training/training/matches_inspection_results/matches_inspection_joined.tsv')
nerpa_results = Path('/home/ilianolhin/git/nerpa2/training/training/nerpa_results/new/temporary')

with approved_matches_table.open('r') as f:
    reader = DictReader(f, delimiter='\t')
    corrected_nrp_ids = [row['NRP variant'].split('#')[0].strip()
                         for row in reader
                         if row['Verdict'] == 'was corrected']

all_matches = []
for result_dir in nerpa_results.iterdir():
    try:
        all_matches.extend(safe_load((result_dir / 'matches_details/matches.yaml').read_text()))
    except FileNotFoundError:
        print(f'No matches for "{result_dir}"')

best_matches = []
for nrp_id in corrected_nrp_ids:
    try:
        best_match = max([match for match in all_matches if match['NRP'] == nrp_id],
                         key=lambda match: match['NormalisedScore'])
    except ValueError:
        print(f'No match for "{nrp_id}"')
        continue
    best_matches.append(best_match)

approved_matches_yaml = Path('/home/ilianolhin/git/nerpa2/training/training/matches_inspection_results/approved_matches.yaml')
approved_matches = safe_load(approved_matches_yaml.read_text())
approved_matches = [match for match in approved_matches if match['NRP'] not in corrected_nrp_ids]
approved_matches.extend(best_matches)

# q: write best_matches to a yaml file
# a:
approved_matches_yaml.write_text(safe_dump(approved_matches))

from prettytable import PrettyTable

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


approved_matches_txt = Path('/home/ilianolhin/git/nerpa2/training/training/matches_inspection_results/approved_matches.txt')
approved_matches_txt.write_text('\n\n'.join(show_match(match) for match in approved_matches))
