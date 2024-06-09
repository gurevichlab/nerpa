#!/usr/bin/env python
import sys
import os
import shutil
import csv

import src.nerpa_pipeline.handle_PCP2 as handle_PCP2
import src.nerpa_pipeline.handle_MT as handle_MT
import src.nerpa_pipeline.handle_E as handle_E
import src.nerpa_pipeline.splitter as splitter
import src.nerpa_pipeline.handle_helper as handle_helper
from collections import defaultdict

from typing import (
    Dict,
    List
)
from src.data_types import (
    BGC_Variant,
    BGC_Module,
    BGC_Module_Modification,
    ResidueScores,
    GeneId,
    dump_bgc_variants
)


# TODO: move the magic constant into the config
MAX_NUM_PARTS = 100


def parse_residue_scores_from_str(scores_line: str) -> ResidueScores:
    residue_scores = {}

    # Split the input string by semicolon to get individual elements
    elements = scores_line.split(';')

    # Iterate through each element and extract the residue and score
    for element in elements:
        parts = element.split('(')
        residue = parts[0].strip()
        score = float(parts[1].replace(')', '').strip())

        # Add the residue and score to the dictionary
        residue_scores[residue] = score

    return residue_scores


def parse_contig_residue_scores(filepath: str) -> Dict[GeneId, List[ResidueScores]]:
    all_residue_scores: Dict[GeneId, List[ResidueScores]] = defaultdict(list)
    with open(filepath, 'r') as rf:
        for line in rf:
            module_name, aa_prediction, scores = line.split('\t')  # TODO: replace scores with proper scores already in this file! I.e. parse aa10, aa34, apply Azat's model, etc.
            ctgorf, domain_id = module_name.rsplit('_', 1)
            # TODO 1: make it just straightforward dump / load in some standard format, e.g. YAML
            # TODO 2: don't use intermediate files at all, parse antiSMASH directly to some Python object
            residue_scores = parse_residue_scores_from_str(scores)
            all_residue_scores[ctgorf].append(residue_scores)
    return all_residue_scores


def build_bgc_assembly_line(orf_modules_names, genome_residue_scores: Dict[GeneId, List[ResidueScores]],
                            dirname) -> List[BGC_Module]:
    '''
    dirname: converted_antiSMASH_v5_outputs
    input_file_name: CONTIGNAME_nrpspredictor2_codes.txt
    orf_part: list of the form ['ctg1_trsI_A1', 'ctg1_trsI_A2']
    '''
    double_orf, double_aa = handle_PCP2.get_double_orfs_and_AA(dirname, orf_modules_names)  # the lists of orfs and modules which can be duplicated
    mt_aa = handle_MT.get_MT_AA(dirname, orf_modules_names)  # the list of modules with methylation
    d_aa = handle_E.get_D_AA(dirname, orf_modules_names)  # the list of modules with epimerization

    bgc_assembly_line: List[BGC_Module] = []
    for ctgorf in orf_modules_names:
        for module_idx, residue_scores in enumerate(genome_residue_scores[ctgorf], start=1):
            module_name = f"{ctgorf}_A{module_idx}"  # FIXME: quite ugly, we expect this particular naming only

            modifications: List[BGC_Module_Modification] = []
            if module_name in mt_aa:
                modifications.append(BGC_Module_Modification.METHYLATION)
            if module_name in d_aa:
                modifications.append(BGC_Module_Modification.EPIMERIZATION)

            cur_bgc_module = BGC_Module(gene_id=ctgorf,
                                        module_idx=int(module_idx),
                                        residue_score=residue_scores,  # TODO: consider renaming in the class definition
                                        modifications=tuple(modifications),  # TODO: consider changing Tuple to List in the class definition
                                        iterative_module=(module_name in double_aa),
                                        iterative_gene=(ctgorf in double_orf))
            bgc_assembly_line.append(cur_bgc_module)
    return bgc_assembly_line


def split_and_reorder(bgc: BGC_Cluster) -> List[BGC_Cluster]:
    return generic_algorithms.list_monad_compose([bgc],
                                                 [splitter.split_by_dist,
                                                  splitter.split_by_single_orf_Starter_TE,
                                                  splitter.split_and_reorder])

def build_bgc_variants(bgc: BGC_Cluster,
                       log: NerpaLogger,
                       residue_scoring_model: Any,
                       config: Any) -> List[BGC_Variant]:  # TODO: replace Any with proper type
    raw_bgc_variants = split_and_reorder(bgc)
    if len(raw_bgc_variants) > config.MAX_VARIANTS_PER_BGC:
        log.info(f'WARNING: Too many parts: {len(parts)}. Keeping first {MAX_NUM_PARTS} of them.')
        raw_bgc_variants = raw_bgc_variants[:MAX_NUM_PARTS]

    iterative_genes, iterative_modules = get_iterative_genes_and_modules(bgc)
                for orf_part in parts[:MAX_NUM_PARTS]:
                    bgc_line = build_bgc_assembly_line(orf_part, genome_residue_scores, dirname)
                    if bgc_line:  # TODO: could it be empty in principle?
                        bgc_variant = BGC_Variant(tentative_assembly_line=bgc_line,
                                                variant_idx=bgc_variant_idx,
                                                genome_id=base_antiSMASHout_name,  # TODO: use proper genome ID from the upstream info
                                                bgc_id=f"bgc#{bgc_variant_idx}")   # TODO: use proper BGC ID from the upstream info (it could be the same for many variants!)
                        bgc_variant_idx += 1
                        all_bgc_variants.append(bgc_variant)

    return all_bgc_variants
