# load approved matches
from csv import DictReader
from yaml import safe_load, safe_dump
from pathlib import Path

approved_matches_table = Path('/home/ilianolhin/git/nerpa2/training/training/matches_inspection_results/matches_inspection_joined.tsv')
all_matches_yaml = Path('/home/ilianolhin/git/nerpa2/training/training/nerpa_results/currently_inspecting/all_matches.yaml')

with approved_matches_table.open('r') as f:
    reader = DictReader(f, delimiter='\t')
    approved_nrp_ids = [row['NRP variant'].split('#')[0].strip()
                        for row in reader
                        if row['Verdict'].startswith('good')]

all_matches = safe_load(all_matches_yaml.read_text())
best_matches = []
for nrp_id in approved_nrp_ids:
    try:
        best_match = max([match for match in all_matches if match['NRP'] == nrp_id],
                         key=lambda match: match['NormalisedScore'])
    except ValueError:
        print(f'No match for "{nrp_id}"')
        continue
    best_matches.append(best_match)

# q: write best_matches to a yaml file
# a:
best_matches_yaml = Path('/home/ilianolhin/git/nerpa2/training/training/matches_inspection_results/approved_matches.yaml')
best_matches_yaml.write_text(safe_dump(best_matches))

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

best_matches_txt = Path('/home/ilianolhin/git/nerpa2/training/training/matches_inspection_results/approved_matches.txt')
best_matches_txt.write_text('\n\n'.join(show_match(match) for match in best_matches))
