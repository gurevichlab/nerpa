from src.matching.matching_types_alignment_step import (
    AlignmentStepType,
    AlignmentStep,
    AlignmentStep_BGC_Module_Info,
    AlignmentStep_NRP_Monomer_Info
)
from src.data_types import NRP_Monomer_Modification
from src.matching.matching_types_alignment import Alignment
from src.matching.matching_types_match import Match
from typing import Iterable, List, Optional, Tuple
from itertools import permutations


def bgc_modules_coincide(module_test: AlignmentStep_BGC_Module_Info,
                         module_approved: AlignmentStep_BGC_Module_Info,
                         indexes_offset: int = 0) -> bool:
    gene_approved = module_approved.gene_id
    if 'ctg' in gene_approved:  # backwards compatibility
        gene_approved = gene_approved.split('_', 1)[1]

    return all([module_test.gene_id == gene_approved,
                module_test.a_domain_idx + indexes_offset == module_approved.a_domain_idx])
                #module_test.modifying_domains == module_approved.modifying_domains,
                #module_test.aa10_code == module_approved.aa10_code,
                #module_test.aa34_code == module_approved.aa34_code])


def filter_out_unknown(modifications: Iterable[NRP_Monomer_Modification]) -> List[NRP_Monomer_Modification]:
    return [mod for mod in modifications if mod != NRP_Monomer_Modification.UNKNOWN]


def nrp_monomers_coincide(monomer_test: AlignmentStep_NRP_Monomer_Info,
                          monomer_approved: AlignmentStep_NRP_Monomer_Info) -> bool:
    return all([monomer_test.chirality == monomer_approved.chirality,
                filter_out_unknown(monomer_test.modifications) == filter_out_unknown(monomer_approved.modifications),
                (monomer_test.rban_name == monomer_approved.rban_name) or
                (monomer_test.rban_name[0] == monomer_approved.rban_name[0] == 'X')])


def steps_coincide(step_test: AlignmentStep, step_approved: AlignmentStep,
                   indexes_offset: int = 0) -> bool:
    result = True
    if step_test.step_type != step_approved.step_type:
        print('Wrong step type')
        result = False

    # compare BGC module info
    if step_test.bgc_module_info is not None and \
            not bgc_modules_coincide(step_test.bgc_module_info, step_approved.bgc_module_info,
                                     indexes_offset):
        print('BGC modules are different')
        result = False

    # compare NRP monomer info
    if step_test.nrp_monomer_info is not None and \
            not nrp_monomers_coincide(step_test.nrp_monomer_info, step_approved.nrp_monomer_info):
        print('NRP monomers are different')
        result = False

    if not result:
        print('Test step:')
        print(step_test)
        print('Approved step:')
        print(step_approved)

    return result


# We compare two alignments by comparing only the MATCH steps
# because the other steps can differ for different versions of Nerpa
def alignments_coincide(alignment1: Alignment,
                        alignment2: Alignment,
                        indexes_offset: int = 0) -> bool:
    min_idx1 = min([step.bgc_module_info.a_domain_idx
                    for step in alignment1
                    if step.bgc_module_info is not None])
    min_idx2 = min([step.bgc_module_info.a_domain_idx
                    for step in alignment2
                    if step.bgc_module_info is not None])
    indexes_offset = min_idx2 - min_idx1
    match_steps1 = [step for step in alignment1 if step.step_type == AlignmentStepType.MATCH]
    match_steps2 = [step for step in alignment2 if step.step_type == AlignmentStepType.MATCH]
    if len(match_steps1) != len(match_steps2):
        print('Different number of MATCH steps')
        return False
    return all(steps_coincide(step1, step2, indexes_offset)
               for step1, step2 in zip(match_steps1,
                                       match_steps2))


def matches_coincide(match1: Match,
                     match2: Match) -> bool:
    joined_alignments1 = [alignment_step
                          for alignment in match1.alignments
                          for alignment_step in alignment]
    joined_alignments2 = [alignment_step
                          for alignment in match2.alignments
                          for alignment_step in alignment]
    return alignments_coincide(joined_alignments1, joined_alignments2)


def find_wrong_match(matches: List[Match], approved_matches: List[Match]) -> Optional[Tuple[Match, Match]]:
    nrp_id_to_match = {match.nrp_variant_info.nrp_id: match
                       for match in matches}

    for approved_match in approved_matches:
        nrp_id = approved_match.nrp_variant_info.nrp_id
        if nrp_id in ['BGC0000296.0',
                      'BGC0001127.0',
                      'BGC0000416.3',
                      'BGC0000307.0',
                      'BGC0001667.2',
                      'BGC0000445.0',
                      'BGC0000339.0',
                      'BGC0000437.0',
                      'BGC0001532.0',
                      'BGC0001667.4',
                      'BGC0002109.2',
                      'BGC0002109.5',
                      'BGC0002109.3',
                      'BGC0002109.1',
                      'BGC0000459.1',
                      'BGC0000450.0']:
            continue
        '''
        if nrp_id in ['BGC0001667.4',
                      'BGC0000296.0',
                      'BGC0001822.1',
                      'BGC0002572.0',
                      'BGC0001127.0',
                      'BGC0001822.3',
                      'BGC0002117.1',
                      'BGC0000416.3',
                      'BGC0000307.0',
                      'BGC0001667.2',
                      'BGC0000383.0',
                      'BGC0002123.11',
                      'BGC0002484.5']:
            continue
        '''
        if nrp_id not in nrp_id_to_match:
            print(f'WARNING: match for {nrp_id} is missing')
            continue
        test_match = nrp_id_to_match[nrp_id]
        if not matches_coincide(test_match, approved_match):
            '''
            with open('/home/ilianolhin/git/nerpa2/training/training/scripts/train_on_groundtruth_matches/wrong_matches.txt', 'a') as f:
                f.write('Wrong match:\n')
                f.write(show_match(match) + '\n\n')
                f.write('Approved match:\n')
                f.write(show_match(approved_match) + '\n\n')
            '''
            return test_match, approved_match
    return None

