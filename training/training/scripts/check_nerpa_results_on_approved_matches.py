from pathlib import Path
from itertools import permutations
import yaml
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


approved_matches_yaml = Path('/home/ilianolhin/git/nerpa2/training/training/matches_inspection_results/approved_matches.yaml')
nerpa_results_dir = Path('/home/ilianolhin/git/nerpa2/training/training/nerpa_results/new/temporary')


def steps_coincide(step1, step2) -> bool:  # step2 is approved
    wrong_field = next((field_name
                        for field_name in ['Alignment_step',
                                           'Modifying_domains',
                                           'NRP_chirality',
                                           'NRP_modifications',
                                           'NRP_residue',
                                           'aa10_code',
                                           'aa34_code']
                        if step1[field_name] != step2[field_name]), None)
    if wrong_field is not None:
        print(wrong_field)
        return False
    if not (step1['rBAN_name'] == step2['rBAN_name'] or step1['rBAN_name'][0] == step2['rBAN_name'][0] == 'X'):
        print('rBAN_name')
        return False
    if not (step1['A-domain_idx'] == step2['A-domain_idx'] == '---' or step1['A-domain_idx'] + 1 == step2['A-domain_idx']):
        print('A-domain_idx')
        return False
    if not (step1['Gene'] == step2['Gene'] == '---' or step1['Gene'] == step2['Gene'].split('_', 1)[1]):
        print('Gene')
        return False
    return all(step1[field_name] == step2[field_name]
               for field_name in ['Alignment_step',
                                  'Modifying_domains',
                                  'NRP_chirality',
                                  'NRP_modifications',
                                  'NRP_residue',
                                  'aa10_code',
                                  'aa34_code'])

def alignments_coincide(alignment1, alignment2) -> bool:
    non_iterating_steps1 = [step for step in alignment1 if step['Alignment_step'] not in ('ITERATE_GENE', 'ITERATE_MODULE')]
    non_iterating_steps2 = [step for step in alignment2 if step['Alignment_step'] not in ('ITERATE_GENE', 'ITERATE_MODULE')]
    return len(non_iterating_steps1) == len(non_iterating_steps2) and\
        all(steps_coincide(step1, step2)
            for step1, step2 in zip(non_iterating_steps1,
                                    non_iterating_steps2))

def matches_coincide(match1, match2) -> bool:
    alignments1 = match1['Alignments']
    alignments2 = match2['Alignments']
    for perm_alignments1 in permutations(alignments1):
        for perm_alignments2 in permutations(alignments2):
            joined_alignments1 = [alignment_step for alignment in perm_alignments1
                                  for alignment_step in alignment]
            joined_alignments2 = [alignment_step for alignment in perm_alignments2
                                  for alignment_step in alignment]
            if alignments_coincide(joined_alignments1, joined_alignments2):
                return True
    return False


nerpa_results_matches = []
for results_dir in nerpa_results_dir.iterdir():
    if not (results_dir / 'matches_details/matches.yaml').exists():
        print(f'{results_dir} does not have matches')
        raise SystemExit
    with open(results_dir / 'matches_details/matches.yaml') as matches_file:
        best_match = max(yaml.safe_load(matches_file),
                         key=lambda match: match['NormalisedScore'],
                         default=None)
        if best_match is None:
            print(f'{results_dir} has no matches')
            raise SystemExit
        nerpa_results_matches.append(best_match)

approved_matches = yaml.safe_load(approved_matches_yaml.read_text())

for i, match in enumerate(nerpa_results_matches):
    nrp_id = match['NRP']
    print(f'{i}. Checking match for {nrp_id}')
    if nrp_id in ['BGC0001098.0',
                  'BGC0002092.0', 'BGC0002092.2', 'BGC0002092.3', 'BGC0002092.5',
                  'BGC0001346.1', 'BGC0001346.2',
                  'BGC0001185.0',
                  'BGC0002564.1','BGC0002564.7',
                  'BGC0002351.0',
                  'BGC0002286.0', 'BGC0002286.1',
                  'BGC0002572.2',
                  'BGC0002370.3',
                  'BGC0002548.0',
                  'BGC0002484.5',
                  'BGC0001035.0',
                  'BGC0002521.0',
                  'BGC0000399.1']:
        print('Skipping for debugging')
        continue
    approved_match = next((approved_match
                           for approved_match in approved_matches
                           if approved_match['NRP'] == nrp_id),
                          None)
    if approved_match is None:
        print(f'NRP {nrp_id} not in approved matches')
        continue
    if not matches_coincide(match, approved_match):
        print(f'Match for {nrp_id} is wrong')
        with open('/home/ilianolhin/git/nerpa2/training/training/scripts/wrong_match.txt', 'w') as f:
            f.write('\n'.join(['Match:',
                               show_match(match),
                               'Approved match:',
                               show_match(approved_match)]))
        break
else:
    print('All matches correct')

