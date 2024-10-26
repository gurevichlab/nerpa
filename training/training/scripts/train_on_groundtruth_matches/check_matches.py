from src.matching.matching_types_alignment_step import (
    AlignmentStepType,
    AlignmentStep,
    AlignmentStep_BGC_Module_Info,
    AlignmentStep_NRP_Monomer_Info
)
from src.matching.matching_types_alignment import Alignment
from src.matching.matching_types_match import Match
from typing import List, Optional, Tuple
from itertools import permutations
from training.training.scripts.train_on_groundtruth_matches.write_results import show_match


def bgc_modules_coincide(module_test: AlignmentStep_BGC_Module_Info,
                         module_approved: AlignmentStep_BGC_Module_Info) -> bool:
    gene_approved = module_approved.gene_id
    if '_' in gene_approved:  # backwards compatibility
        gene_approved = gene_approved.split('_', 1)[1]

    return all([module_test.gene_id == gene_approved,
                module_test.a_domain_idx == module_approved.a_domain_idx,
                module_test.modifying_domains == module_approved.modifying_domains,
                module_test.aa10_code == module_approved.aa10_code,
                module_test.aa34_code == module_approved.aa34_code])


def nrp_monomers_coincide(monomer_test: AlignmentStep_NRP_Monomer_Info,
                          monomer_approved: AlignmentStep_NRP_Monomer_Info) -> bool:
    return all([monomer_test.chirality == monomer_approved.chirality,
                monomer_test.modifications == monomer_approved.modifications,
                monomer_test.rban_name == monomer_approved.rban_name])


def steps_coincide(step_test: AlignmentStep, step_approved: AlignmentStep) -> bool:
    if step_test.step_type != step_approved.step_type:
        print('Wrong step type')
        return False

    # compare BGC module info
    if step_test.bgc_module_info is not None and \
            not bgc_modules_coincide(step_test.bgc_module_info, step_approved.bgc_module_info):
        print('BGC module')
        return False

    # compare NRP monomer info
    if step_test.nrp_monomer_info is not None and \
            not nrp_monomers_coincide(step_test.nrp_monomer_info, step_approved.nrp_monomer_info):
        print('NRP monomer')
        return False

    return True


# We compare two alignments by comparing only the MATCH steps
# because the other steps can differ for different versions of Nerpa
def alignments_coincide(alignment1: Alignment,
                        alignment2: Alignment) -> bool:
    match_steps1 = [step for step in alignment1 if step.step_type == AlignmentStepType.MATCH]
    match_steps2 = [step for step in alignment2 if step.step_type == AlignmentStepType.MATCH]
    return len(match_steps1) == len(match_steps2) and \
        all(steps_coincide(step1, step2)
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

