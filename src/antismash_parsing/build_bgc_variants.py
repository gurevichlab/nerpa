from typing import (
    Any,
    Callable,
    Dict,
    List
)
from src.antismash_parsing.antismash_parser_types import (
    A_Domain,
    BGC_Cluster,
    DomainType,
    Gene,
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
    BGC_Fragment,
    LogProb
)
from src.monomer_names_helper import MonomerNamesHelper, antiSMASH_MonomerName, MonomerResidue
from src.pipeline.logger import NerpaLogger
from src.config import antiSMASH_Parsing_Config
from src.antismash_parsing.bgcs_split_and_reorder import split_and_reorder
from src.generic.string import hamming_distance
from itertools import chain
import pandas as pd
from collections import defaultdict


def get_residue_scores(a_domain: A_Domain,
                       residue_scoring_model: Callable[[pd.DataFrame], Dict[MonomerResidue, LogProb]],
                       monomer_names_helper: MonomerNamesHelper,
                       config: antiSMASH_Parsing_Config) -> Dict[MonomerResidue, LogProb]:
    def svm_score(aa_name: MonomerResidue, svm_level_prediction: SVM_Prediction) -> float:
        svm_prediction_substrates = {monomer_names_helper.parsed_name(monomer_name, 'antismash').residue
                                     for monomer_name in svm_level_prediction.substrates}
        if aa_name not in config.SVM_SUBSTRATES:
            return config.SVM_NOT_SUPPORTED_SCORE
        if aa_name in svm_prediction_substrates:
            return svm_level_prediction.score
        else:
            return config.SVM_NO_PREDICTION_SCORE

    def similarity_score(aa_code1: str, aa_code2: str) -> float:
        return 1.0 - hamming_distance(aa_code1, aa_code2) / len(aa_code1)

    scoring_table = pd.DataFrame([], columns=config.SCORING_TABLE_COLUMNS).set_index(config.SCORING_TABLE_INDEX)
    for aa_name, aa10_codes in config.KNOWN_AA10_CODES.items():
        aa34_codes = config.KNOWN_AA34_CODES[aa_name]
        aa_code_scores = [max(similarity_score(aa_code, known_aa_code)
                              for known_aa_code in known_aa_codes)
                          for aa_code, known_aa_codes in [(a_domain.aa10, aa10_codes),
                                                          (a_domain.aa34, aa34_codes)]]
        svm_scores = [svm_score(aa_name, a_domain.svm[level])
                      for level in SVM_LEVEL]
        scoring_table.loc[aa_name] = aa_code_scores + svm_scores

    return residue_scoring_model(scoring_table)


def build_gene_assembly_line(gene: Gene,
                             residue_scoring_model: Any,
                             monomer_names_helper: MonomerNamesHelper,
                             config: antiSMASH_Parsing_Config) -> List[BGC_Module]:
    iterative_gene = is_iterative_gene(gene)
    iterative_modules_idxs = get_iterative_modules_idxs(gene)
    built_modules = []
    a_domain_idx = 0  # indexes start from 1 for backward compatibility. Maybe change this in the future
    for module_idx, module in enumerate(gene.modules):
        if module.a_domain is None:
            continue
        a_domain_idx += 1

        residue_scores = get_residue_scores(module.a_domain,
                                            residue_scoring_model,
                                            monomer_names_helper,
                                            config)

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
                            monomer_names_helper: MonomerNamesHelper,
                            config: antiSMASH_Parsing_Config) -> List[BGC_Module]:
    return list(chain.from_iterable(build_gene_assembly_line(gene,
                                                             residue_scoring_model,
                                                             monomer_names_helper,
                                                             config)
                                    for gene in bgc_genes))


def build_bgc_fragments(raw_fragmented_bgc: List[BGC_Cluster],
                        residue_scoring_model: Any,
                        monomer_names_helper: MonomerNamesHelper,
                        config: antiSMASH_Parsing_Config) -> List[BGC_Fragment]:
    return [build_bgc_assembly_line(bgc_cluster.genes,
                                    residue_scoring_model,
                                    monomer_names_helper,
                                    config)
            for bgc_cluster in raw_fragmented_bgc]


def build_bgc_variants(bgc: BGC_Cluster,
                       residue_scoring_model: Any,
                       monomer_names_helper: MonomerNamesHelper,
                       config: antiSMASH_Parsing_Config,
                       log: NerpaLogger) -> List[BGC_Variant]:  # TODO: replace Any with proper type
    raw_fragmented_bgcs = split_and_reorder(bgc, config)
    if len(raw_fragmented_bgcs) > config.MAX_VARIANTS_PER_BGC:
        log.info(f'WARNING: Too many parts: {len(raw_fragmented_bgcs)}. Keeping first {config.MAX_VARIANTS_PER_BGC} of them.')
        raw_fragmented_bgcs = raw_fragmented_bgcs[:config.MAX_VARIANTS_PER_BGC]


    return [BGC_Variant(genome_id=bgc.genome_id,
                        variant_idx=idx,
                        bgc_idx=bgc.bgc_idx,
                        fragments=build_bgc_fragments(raw_fragmented_bgc,
                                                      residue_scoring_model,
                                                      monomer_names_helper,
                                                      config),
                        has_pks_domains=bgc.has_pks_domains())
            for idx, raw_fragmented_bgc in enumerate(raw_fragmented_bgcs)]
