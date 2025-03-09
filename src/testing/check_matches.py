from src.matching.alignment_step_type import (
    AlignmentStep,
    AlignmentStep_BGC_Module_Info
)
from src.data_types import NRP_Monomer_Modification
from src.matching.alignment_type import Alignment
from src.matching.match_type import Match
from src.matching.hmm_auxiliary_types import DetailedHMMEdgeType
from src.rban_parsing.rban_monomer import rBAN_Monomer
from typing import Iterable, List, Optional, Tuple, NamedTuple


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


def nrp_monomers_coincide(monomer_test: rBAN_Monomer,
                          monomer_approved: rBAN_Monomer) -> bool:
    return (monomer_test.rban_name == monomer_approved.rban_name and
            monomer_test.chirality == monomer_approved.chirality) or \
        (monomer_test.rban_name[0] == monomer_approved.rban_name[0] == 'X')


def steps_coincide(step_test: AlignmentStep, step_approved: AlignmentStep,
                   indexes_offset: int = 0) -> bool:
    result = True
    if step_test.step_type != step_approved.step_type:
        print('Wrong step type')
        result = False

    # compare BGC module info
    if step_test.bgc_module is not None and \
            not bgc_modules_coincide(step_test.bgc_module, step_approved.bgc_module,
                                     indexes_offset):
        print('BGC modules are different')
        result = False

    # compare NRP monomer info
    if step_test.nrp_monomer is not None and \
            not nrp_monomers_coincide(step_test.nrp_monomer, step_approved.nrp_monomer):
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
    min_idx1 = min([step.bgc_module.a_domain_idx
                    for step in alignment1
                    if step.bgc_module is not None])
    min_idx2 = min([step.bgc_module.a_domain_idx
                    for step in alignment2
                    if step.bgc_module is not None])
    indexes_offset = min_idx2 - min_idx1
    match_steps1 = [step for step in alignment1 if step.step_type == DetailedHMMEdgeType.MATCH]
    match_steps2 = [step for step in alignment2 if step.step_type == DetailedHMMEdgeType.MATCH]
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


class MatchPair(NamedTuple):
    test_match: Match
    approved_match: Match

class CheckResults(NamedTuple):
    missing_cnt: int
    wrong_matches: List[MatchPair]


def find_wrong_matches(matches: List[Match], approved_matches: List[Match]) -> CheckResults:
    missing_cnt = 0
    wrong_matches = []
    for approved_match in approved_matches:
        nrp_id = approved_match.nrp_variant_id.nrp_id
        bgc_id = nrp_id.split('.')[0]
        try:
            test_match = next(match for match in matches
                              if match.nrp_variant_id.nrp_id == nrp_id
                              and match.bgc_variant_id.bgc_id.genome_id == bgc_id)
        except StopIteration:
            print(f'WARNING: match for {nrp_id} is missing')
            missing_cnt += 1
            continue
        print('Testing ', nrp_id)
        if not matches_coincide(test_match, approved_match):
            wrong_matches.append(MatchPair(test_match, approved_match))

    return CheckResults(missing_cnt, wrong_matches)
