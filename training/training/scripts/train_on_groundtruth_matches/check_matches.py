from typing import List, Optional, Tuple
from itertools import permutations
from training.training.scripts.train_on_groundtruth_matches.write_results import show_match


def steps_coincide(step_test, step_approved) -> bool:  # step2 is approved
    wrong_field = next((field_name
                        for field_name in ['Alignment_step',
                                           'Modifying_domains',
                                           'NRP_chirality',
                                           #'NRP_modifications',
                                           #'NRP_residue',
                                           'aa10_code',
                                           'aa34_code']
                        if step_test[field_name] != step_approved[field_name]), None)
    if wrong_field is not None:
        print(wrong_field)
        return False
    if not (step_test['rBAN_name'] == step_approved['rBAN_name'] or step_test['rBAN_name'][0] == step_approved['rBAN_name'][0] == 'X'):
        print('rBAN_name')
        return False
    if not (step_test['A-domain_idx'] == step_approved['A-domain_idx'] == '---' or step_test['A-domain_idx'] == step_approved['A-domain_idx']):
        print('A-domain_idx')
        return False
    if not (step_test['Gene'] == step_approved['Gene'] or ('_' in step_approved['Gene'] and step_test['Gene'] == step_approved['Gene'].split('_', 1)[1])):
        print('Gene')
        return False
    return all(step_test[field_name] == step_approved[field_name]
               for field_name in ['Alignment_step',
                                  'Modifying_domains',
                                  'NRP_chirality',
                                  #'NRP_modifications',
                                  'aa10_code',
                                  'aa34_code'])


def alignments_coincide(alignment1, alignment2) -> bool:
    match_steps1 = [step for step in alignment1 if step['Alignment_step'] == 'MATCH']
    match_steps2 = [step for step in alignment2 if step['Alignment_step'] == 'MATCH']
    return len(match_steps1) == len(match_steps2) and \
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


def find_wrong_match(matches: List[dict], approved_matches: List[dict]) -> Optional[Tuple[str, str]]:
    nrp_id_to_match = {match['NRP']: match for match in matches}
    nrp_id_to_approved_match = {match['NRP']: match for match in approved_matches}

    for match in matches:
        if match['NRP'] not in nrp_id_to_approved_match:
            continue
        approved_match = nrp_id_to_approved_match[match['NRP']]
        if not matches_coincide(match, approved_match):
            #return match['NRP'], show_match(match)

            with open('/home/ilianolhin/git/nerpa2/training/training/scripts/train_on_groundtruth_matches/wrong_matches.txt', 'a') as f:
                f.write('Wrong match:\n')
                f.write(show_match(match) + '\n\n')
                f.write('Approved match:\n')
                f.write(show_match(approved_match) + '\n\n')

