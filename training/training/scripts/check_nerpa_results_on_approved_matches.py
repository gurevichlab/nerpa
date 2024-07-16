from pathlib import Path
from itertools import permutations
import yaml


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
                                           'aa34_code',
                                           'rBAN_idx',
                                           'rBAN_name']
                        if step1[field_name] != step2[field_name]), None)
    if wrong_field:
        print(wrong_field)
        return False
    if step1['A-domain_idx'] + 1 != step2['A-domain_idx']:
        print('A-domain_idx')
        return False
    if step1['Gene'] != step2['Gene'].split('_', 1)[1]:
        print('Gene')
        return False
    return all(step1[field_name] == step2[field_name]
               for field_name in ['Alignment_step',
                                  'Modifying_domains',
                                  'NRP_chirality',
                                  'NRP_modifications',
                                  'NRP_residue',
                                  'aa10_code',
                                  'aa34_code',
                                  'rBAN_idx',
                                  'rBAN_name'])

def alignments_coincide(alignment1, alignment2) -> bool:
    return len(alignment1) == len(alignment2) and all(steps_coincide(step1, step2)
                                                      for step1, step2 in zip(alignment1, alignment2))

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
    with open(results_dir / 'matches_details/matches.yaml') as matches_file:
        best_match = max(yaml.safe_load(matches_file),
                         key=lambda match: match['NormalisedScore'])
        nerpa_results_matches.append(best_match)

approved_matches = yaml.safe_load(approved_matches_yaml.read_text())

for match in nerpa_results_matches:
    nrp_id = match['NRP']
    approved_match = next((approved_match
                           for approved_match in approved_matches
                           if approved_match['NRP'] == nrp_id),
                          None)
    if approved_match is None:
        print(f'NRP {nrp_id} not in approved matches')
        break
    if not matches_coincide(match, approved_match):
        print(f'Match for {nrp_id} is wrong')
        break
else:
    print('All matches correct')

