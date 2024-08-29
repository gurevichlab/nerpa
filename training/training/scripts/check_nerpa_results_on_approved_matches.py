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
                                           #'NRP_modifications',
                                           #'NRP_residue',
                                           'aa10_code',
                                           'aa34_code']
                        if step1[field_name] != step2[field_name]), None)
    if wrong_field is not None:
        print(wrong_field)
        return False
    if not (step1['rBAN_name'] == step2['rBAN_name'] or step1['rBAN_name'][0] == step2['rBAN_name'][0] == 'X'):
        print('rBAN_name')
        return False
    if not (step1['A-domain_idx'] == step2['A-domain_idx'] == '---' or step1['A-domain_idx'] == step2['A-domain_idx']):
        print('A-domain_idx')
        return False
    if not (step1['Gene'] == step2['Gene'] or ('_' in step2['Gene'] and step1['Gene'] == step2['Gene'].split('_', 1)[1])):
        print('Gene')
        return False
    return all(step1[field_name] == step2[field_name]
               for field_name in ['Alignment_step',
                                  'Modifying_domains',
                                  'NRP_chirality',
                                  #'NRP_modifications',
                                  'aa10_code',
                                  'aa34_code'])

def alignments_coincide(alignment1, alignment2) -> bool:
    match_steps1 = [step for step in alignment1 if step['Alignment_step'] == 'MATCH']
    match_steps2 = [step for step in alignment2 if step['Alignment_step'] == 'MATCH']
    return len(match_steps1) == len(match_steps2) and\
        all(steps_coincide(step1, step2)
            for step1, step2 in zip(match_steps1,
                                    match_steps2))

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

approved_matches = yaml.safe_load(approved_matches_yaml.read_text())
print(f'{len(approved_matches)} approved matches in total')

approved_nrp_ids = {match['NRP'] for match in approved_matches}

nerpa_results_matches = []
for results_dir in nerpa_results_dir.iterdir():
    if not results_dir.name in approved_nrp_ids:
        continue
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

wrong_matches = []
for i, match in enumerate(nerpa_results_matches):
    nrp_id = match['NRP']
    print(f'{i}. Checking match for {nrp_id}')
    approved_match = next((approved_match
                           for approved_match in approved_matches
                           if approved_match['NRP'] == nrp_id),
                          None)
    if approved_match is None:
        print(f'NRP {nrp_id} not in approved matches')
        continue
    '''
    if nrp_id in ['BGC0002109.4',
                  'BGC0001421.0',
                  'BGC0002188.0',
                  'BGC0000374.4',
                  'BGC0000379.0',
                  'BGC0001955.2',
                  'BGC0000439.0',
                  'BGC0001214.2',
                  'BGC0001230.7',
                  'BGC0001955.3',
                  'BGC0001214.4',
                  ]:
    
        print('Skipping for debugging')
        continue
    '''
    if nrp_id in ['BGC0000374.4',  # corrected
                  'BGC0000379.0',  # to check
                  'BGC0000439.0',  # to check
                  'BGC0001230.7',  # to check
                  'BGC0002415.2',  # corrected
                  'BGC0001214.4',  # corrected
                  'BGC0002415.1'  # corrected
                  ]:
        print('Skipping for debugging')
        continue
    if not matches_coincide(match, approved_match):
        print(f'Match for {nrp_id} is wrong')
        wrong_matches.append(nrp_id)
        with open('/home/ilianolhin/git/nerpa2/training/training/scripts/wrong_match.txt', 'w') as f:
            f.write('\n'.join(['Match:',
                               show_match(match),
                               'Approved match:',
                               show_match(approved_match),
                               '\n\n']))
        break
else:
    print('All matches checked')
    if not wrong_matches:
        print('Everything correct')
    else:
        print('Wrong matches:')
        print(wrong_matches)
        print('Check /home/ilianolhin/git/nerpa2/training/training/scripts/wrong_match.txt')


