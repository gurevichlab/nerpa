from typing import List, Tuple, Optional, TypeVar, Callable, Hashable, Iterable

from src.antismash_parsing.bgc_variant_types import BGC_Variant, A_Domain_ID
from src.matching.alignment_type import Alignment
from src.monomer_names_helper import MonomerNamesHelper
from src.matching.match_type import Match
from src.generic.combinatorics import is_subsequence, filter_unique, remove_runs_of_equal_elements, \
    split_into_iterations
from itertools import pairwise, groupby, chain, takewhile
from more_itertools import peekable

from src.testing.testing_types import TestMatch
from src.testing.simplified_alignment import split_into_bgc_iterations, SimplifiedAlignment


def remove_gene_module_iterations(a_domains: List[A_Domain_ID]) -> List[A_Domain_ID]:
    result = []
    for gene_id, a_domains_group in groupby(a_domains, lambda aid: aid.gene_id):
        gene_iterations = split_into_iterations(a_domains_group)
        result.extend(max(gene_iterations, key=len))

    # gene iterations removed, now remove module iterations
    return list(remove_runs_of_equal_elements(result))

def bgc_variant_match_compatible(bgc_variant: BGC_Variant,
                                 match: TestMatch) -> bool:
    def aids_from_alignment(alignment: SimplifiedAlignment) -> List[A_Domain_ID]:
        return [a_domain_id
                for a_domain_id, _ in alignment
                if a_domain_id is not None]

    bgc_a_domains = [A_Domain_ID.from_module(module)
                     for module in bgc_variant.modules]

    alignment_a_domain_lists = [aids_from_alignment(bgc_iteration)
                                for bgc_iteration in split_into_bgc_iterations(match.true_alignment)]

    # remove runs of equal modules from alignment_modules (those are modules/genes iterations)
    alignment_a_domain_lists = [remove_gene_module_iterations(alignment_a_domain_list)
                                for alignment_a_domain_list in alignment_a_domain_lists]

    if not all(is_subsequence(alignment_a_domain_list, bgc_a_domains)
               for alignment_a_domain_list in alignment_a_domain_lists):
        return False
    return True


# def fix_bgc_variant(bgc_variant: BGC_Variant) -> BGC_Variant:
#     if (min_idx := min(module.fragment_idx for module in bgc_variant.modules)) > 0:
#         for module in bgc_variant.modules:
#             module.fragment_idx -= min_idx
#
#     return bgc_variant

# def fix_indexing_in_alignment(alignment: Alignment):
#     if min(step.bgc_module.a_domain_idx
#            for step in alignment
#            if step.bgc_module is not None) == 1:
#         for step in alignment:
#             if step.bgc_module is not None:
#                 step.bgc_module = step.bgc_module._replace(a_domain_idx=step.bgc_module.a_domain_idx - 1)
#
#
# def fix_gene_names_in_alignment(alignment: Alignment):
#     if all(step.bgc_module.gene_id.startswith('ctg')
#            for step in alignment
#            if step.bgc_module is not None) == 1:
#         for step in alignment:
#             if step.bgc_module is not None:
#                 step.bgc_module = step.bgc_module._replace(gene_id=GeneId(step.bgc_module.gene_id.split('_')[1]))
#
#
# def fix_residue_names_in_alignment(alignment: Alignment,
#                                    monomer_names_helper: MonomerNamesHelper):
#     for step in alignment:
#         if step.nrp_monomer is not None:
#             parsed_monomer = monomer_names_helper.parsed_name(step.nrp_monomer.rban_name, name_format='norine')
#             step.nrp_monomer = step.nrp_monomer._replace(residue=parsed_monomer.residue)
#
#
# def fix_old_style_multiple_alignments(match: Match) -> Match:
#     # Previously, different parts of BGC could correspond to multiple alignments
#     # Here, I join all alignments which don't have common modules into one alignment
#     if len(match.alignments) == 1:
#         return match
#     new_alignments = []
#     new_alignment = []
#     for alignment in match.alignments:
#         prev_modules = [(step.bgc_module.gene_id, step.bgc_module.a_domain_idx)
#                         for step in new_alignment
#                         if step.bgc_module is not None]
#         alignment_modules = [(step.bgc_module.gene_id, step.bgc_module.a_domain_idx)
#                              for step in alignment
#                              if step.bgc_module is not None]
#         if not set(prev_modules) & set(alignment_modules):
#             new_alignment.extend(alignment)
#         else:
#             new_alignments.append(new_alignment)
#             new_alignment = alignment
#     new_alignments.append(new_alignment)
#     match.alignments = new_alignments
#     return match
#
#
# def fix_alignment(alignment: Alignment,
#                   monomer_names_helper: MonomerNamesHelper) -> Alignment:
#     fix_indexing_in_alignment(alignment)
#     fix_gene_names_in_alignment(alignment)
#     fix_residue_names_in_alignment(alignment,
#                                    monomer_names_helper)
#     return alignment
#
#
# def fix_match(match: Match,
#               monomer_names_helper: MonomerNamesHelper) -> Match:
#     for alignment in match.alignments:
#         fix_alignment(alignment, monomer_names_helper)
#     fix_old_style_multiple_alignments(match)
#     return match
#
