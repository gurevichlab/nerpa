from typing import (
    Any,
    List
)
from src.antismash_parsing.antismash_parser_types import (
    A_Domain,
    BGC_Cluster,
    DomainType,
    Gene,
    MonomerResidue,
    SVM_LEVEL,
    SVM_Prediction
)
from src.antismash_parsing.determine_modifications import (
    get_iterative_modules_idxs,
    is_iterative_gene
)
from src.data_types import (
    BGC_Variant,
    BGC_Module,
    BGC_Module_Modification,
    ResidueScores,
)
from src.pipeline.logger import NerpaLogger
from src.config import antiSMASH_Parsing_Config
from src.antismash_parsing.bgcs_split_and_reorder import split_and_reorder
from src.generic.string import hamming_distance
from itertools import chain
import pandas as pd

#TODO: needs typing for residue_scoring_model and config
def get_residue_scores(a_domain: A_Domain,
                       residue_scoring_model: Any,
                       config: antiSMASH_Parsing_Config) -> ResidueScores:
    def svm_score(residue: MonomerResidue, svm_level_prediction: SVM_Prediction) -> float:
        if residue not in config.SVM_SUBSTRATES:
            return -1.0
        if residue in svm_level_prediction.monomer_residues:
            return svm_level_prediction.score
        else:
            return 0.0

    def similarity_score(aa_code1: str, aa_code2: str) -> float:
        return 1.0 - hamming_distance(aa_code1, aa_code2) / len(aa_code1)

    scoring_table = pd.DataFrame([], columns=config.SCORING_TABLE_COLUMNS).set_index(config.SCORING_TABLE_INDEX)
    for residue, aa10_codes in config.KNOWN_AA10_CODES.items():
        aa_34_codes = config.KNOWN_AA34_CODES[residue]
        aa_code_scores = [max(similarity_score(aa_code, known_aa_code)
                              for known_aa_code in known_aa_codes)
                          for aa_code, known_aa_codes in [(a_domain.aa10, aa10_codes),
                                                          (a_domain.aa34, aa_34_codes)]]
        svm_scores = [svm_score(residue, a_domain.svm[level])
                      for level in SVM_LEVEL]
        scoring_table.loc[residue] = aa_code_scores + svm_scores

    return residue_scoring_model(scoring_table)


def build_gene_assembly_line(gene: Gene,
                             residue_scoring_model: Any,
                             config: antiSMASH_Parsing_Config) -> List[BGC_Module]:
    iterative_gene = is_iterative_gene(gene)
    iterative_modules_idxs = get_iterative_modules_idxs(gene)
    built_modules = []
    a_domain_idx = -1
    for module_idx, module in enumerate(gene.modules):
        if module.a_domain is None:
            continue
        a_domain_idx += 1

        residue_scores = get_residue_scores(module.a_domain, residue_scoring_model, config)

        modifications = []
        if DomainType.MT in module.domains_sequence:
            modifications.append(BGC_Module_Modification.METHYLATION)
        if DomainType.E in module.domains_sequence or \
                (module_idx < len(gene.modules)-1 and DomainType.C_DUAL in gene.modules[module_idx+1].domains_sequence):
            modifications.append(BGC_Module_Modification.EPIMERIZATION)

        built_modules.append(BGC_Module(gene_id=gene.gene_id,
                                        a_domain_idx=a_domain_idx,
                                        residue_score=residue_scores,
                                        modifications=tuple(modifications),
                                        iterative_module=module_idx in iterative_modules_idxs,
                                        iterative_gene=iterative_gene,
                                        aa10_code=module.a_domain.aa10,
                                        aa34_code=module.a_domain.aa34))

    return built_modules


def build_bgc_assembly_line(bgc_genes: List[Gene],
                            residue_scoring_model: Any,
                            config: antiSMASH_Parsing_Config) -> List[BGC_Module]:
    return list(chain.from_iterable(build_gene_assembly_line(gene,
                                                             residue_scoring_model,
                                                             config)
                                    for gene in bgc_genes))


def build_bgc_variants(bgc: BGC_Cluster,
                       log: NerpaLogger,
                       residue_scoring_model: Any,
                       config: antiSMASH_Parsing_Config) -> List[BGC_Variant]:  # TODO: replace Any with proper type
    raw_bgc_variants = split_and_reorder(bgc, config)
    if len(raw_bgc_variants) > config.MAX_VARIANTS_PER_BGC:
        log.info(f'WARNING: Too many parts: {len(raw_bgc_variants)}. Keeping first {config.MAX_VARIANTS_PER_BGC} of them.')
        raw_bgc_variants = raw_bgc_variants[:config.MAX_VARIANTS_PER_BGC]

    return [BGC_Variant(genome_id=bgc.genome_id,
                        variant_idx=idx,
                        bgc_idx=bgc.bgc_idx,
                        tentative_assembly_line=build_bgc_assembly_line(raw_bgc_variant.genes,
                                                                        residue_scoring_model,
                                                                        config))
            for idx, raw_bgc_variant in enumerate(raw_bgc_variants)]
