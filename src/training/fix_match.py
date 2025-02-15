from typing import List
from src.matching.matching_types_alignment import Alignment
from src.data_types import GeneId
from src.monomer_names_helper import monomer_names_helper
from src.matching.matching_types_match import Match
from src.data_types import GeneId, BGC_Variant
from src.generic.combinatorics import is_subsequence, filter_unique
from itertools import pairwise

def fix_indexing_in_alignment(alignment: Alignment):
    if min(step.bgc_module.a_domain_idx
           for step in alignment
           if step.bgc_module is not None) == 1:
        for step in alignment:
            if step.bgc_module is not None:
                step.bgc_module = step.bgc_module._replace(a_domain_idx=step.bgc_module.a_domain_idx - 1)

def fix_gene_names_in_alignment(alignment: Alignment):
    if all(step.bgc_module.gene_id.startswith('ctg')
           for step in alignment
           if step.bgc_module is not None) == 1:
        for step in alignment:
            if step.bgc_module is not None:
                step.bgc_module = step.bgc_module._replace(gene_id=GeneId(step.bgc_module.gene_id.split('_')[1]))

def fix_residue_names_in_alignment(alignment: Alignment):
    for step in alignment:
        if step.nrp_monomer is not None:
            parsed_monomer = monomer_names_helper.parsed_name(step.nrp_monomer.rban_name, name_format='norine')
            step.nrp_monomer = step.nrp_monomer._replace(residue=parsed_monomer.residue)


def fix_old_style_multiple_alignments(match: Match) -> Match:
    # Previously, different parts of BGC could correspond to multiple alignments
    # Here, I join all alignments which don't have common modules into one alignment
    if len(match.alignments) == 1:
        return match
    new_alignments = []
    new_alignment = []
    for alignment in match.alignments:
        prev_modules = [(step.bgc_module.gene_id, step.bgc_module.a_domain_idx)
                        for step in new_alignment
                        if step.bgc_module is not None]
        alignment_modules = [(step.bgc_module.gene_id, step.bgc_module.a_domain_idx)
                      for step in alignment
                      if step.bgc_module is not None]
        if not set(prev_modules) & set(alignment_modules):
            new_alignment.extend(alignment)
        else:
            new_alignments.append(new_alignment)
            new_alignment = alignment
    new_alignments.append(new_alignment)
    match.alignments = new_alignments
    return match


def fix_alignment(alignment: Alignment) -> Alignment:
    fix_indexing_in_alignment(alignment)
    fix_gene_names_in_alignment(alignment)
    fix_residue_names_in_alignment(alignment)
    return alignment


def fix_match(match: Match) -> Match:
    for alignment in match.alignments:
        fix_alignment(alignment)
    fix_old_style_multiple_alignments(match)
    return match


def bgc_variant_match_compatible(bgc_variant: BGC_Variant,
                                 match: Match) -> bool:
    bgc_modules = [(module.gene_id, module.a_domain_idx)
                   for module in bgc_variant.modules]
    for alignment in match.alignments:
        alignment_modules = [(step.bgc_module.gene_id, step.bgc_module.a_domain_idx)
                             for step in alignment
                             if step.bgc_module is not None
                             and step.nrp_monomer is not None]
        # remove runs of equal modules from alignment_modules (those are modules/genes iterations)
        alignment_modules = list(filter_unique(alignment_modules))

        if not is_subsequence(alignment_modules, bgc_modules):
            return False
    return True


def fix_bgc_variant(bgc_variant: BGC_Variant) -> BGC_Variant:
    if (min_idx := min(module.fragment_idx for module in bgc_variant.modules)) > 0:
        for module in bgc_variant.modules:
            module.fragment_idx -= min_idx

    return bgc_variant