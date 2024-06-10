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

#TODO: needs typing for residue_scoring_model and config
def get_residue_scores(a_domain: A_Domain,
                       known_aa10_codes: Dict[MonomerResidue, List[str]],
                       known_aa34_codes: Dict[MonomerResidue, List[str]],
                       residue_scoring_model: Any,
                       config: Any) -> ResidueScores:
    def svm_score(residue: MonomerResidue, svm_level_prediction: SVM_Prediction) -> float:
        if residue not in config.SVM_SUBSTRATES:
            return -1.0
        if residue in svm_level_prediction.monomer_residues:
            return svm_level_prediction.score
        else:
            return 0.0

    scoring_table = pd.DataFrame([], columns=config.SCORING_TABLE_COLUMNS).set_index(config.SCORING_TABLE_INDEX)
    for residue, aa10_codes in known_aa10_codes.items():
        aa_34_codes = known_aa34_codes[residue]
        aa_code_scores = [min_hamming_distance(aa10_code, known_aa10_codes),
                          min_hamming_distance(aa34_code, known_aa34_codes)]
        svm_scores = [svm_score(residue, a_domain.svm[level])
                      for level in SVM_LEVEL]
        scoring_table.loc[residue] = aa_code_scores + svm_scores

    return residue_scoring_model(scoring_table)


def build_gene_assembly_line(gene: Gene,
                             residue_scoring_model: Any) -> List[BGC_Module]:
    iterative_gene = is_iterative_gene(gene)
    iterative_modules_idxs = get_iterative_modules_idxs(gene)
    built_modules = []
    for module_idx, module in enumerate(gene.modules):
        if module.a_domain is None:
            continue

        residue_scores = get_residue_scores(module.a_domain, residue_scoring_model)

        modifications = []
        if DomainType.MT in module.domains_sequence:
            modifications.append(BGC_Module_Modification.METHYLATION)
        if DomainType.E in module.domains_sequence or DomainType.C_DUAL in module.domains_sequence:
            modifications.append(BGC_Module_Modification.EPIMERIZATION)

        built_modules.append(BGC_Module(gene_id=gene.gene_id,
                                        module_idx=module_idx,
                                        residue_score=residue_scores,
                                        modifications=tuple(modifications),
                                        iterative_module=module_idx in iterative_modules_idxs,
                                        iterative_gene=iterative_gene))

    return built_modules


def build_bgc_assembly_line(bgc_genes: List[Gene],
                            residue_scoring_model: Any) -> List[BGC_Module]:
    return list(chain(*build_gene_assembly_line(gene, residue_scoring_model)
                      for gene in bgc_genes))

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
        log.info(f'WARNING: Too many parts: {len(raw_bgc_variants)}. Keeping first {MAX_NUM_PARTS} of them.')
        raw_bgc_variants = raw_bgc_variants[:MAX_NUM_PARTS]

    return [BGC_Variant(tentative_assembly_line=build_bgc_assembly_line(raw_bgc_variant, residue_scoring_model),
                        variant_idx=idx,
                        genome_id=bgc.genome_id,
                        bgc_id=bgc.bgc_idx)
            for idx, raw_bgc_variant in enumerate(raw_bgc_variants)]
