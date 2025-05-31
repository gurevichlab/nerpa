from copy import deepcopy
from math import log
from typing import (
    Any,
    Callable,
    Dict,
    List,
    NamedTuple,
    Optional
)
from src.aa_specificity_prediction_model.model_wrapper import ModelWrapper
from src.aa_specificity_prediction_model.specificity_prediction_helper import SpecificityPredictionHelper
from src.antismash_parsing.antismash_parser_types import (
    A_Domain,
    BGC_Cluster,
    DomainType,
    Gene,
    SVM_LEVEL,
    SVM_Prediction, BGC_Module_ID, Fragmented_BGC_Cluster
)
from src.antismash_parsing.determine_modifications import (
    get_modules_modifications,
    get_iterative_genes,
    get_iterative_modules_ids
)
from src.data_types import (
    AA34,
    BGC_ID,
    BGC_Variant_ID,
    BGC_Variant,
    BGC_Module,
    LogProb
     )
from src.monomer_names_helper import MonomerNamesHelper, MonomerResidue
from src.pipeline.logger import NerpaLogger
from src.config import antiSMASH_Processing_Config, SpecificityPredictionConfig
from src.antismash_parsing.bgcs_split_and_reorder import split_and_reorder
from src.generic.string import hamming_distance
from src.generic.functional import cached_by_key
from src.antismash_parsing.genomic_context import get_modules_genomic_context
from itertools import chain


def get_residue_scores(a_domain: A_Domain,
                       specificity_prediction_helper: SpecificityPredictionHelper) -> Dict[MonomerResidue, LogProb]:
    #apriori_prob = specificity_prediction_helper.config.APRIORI_RESIDUE_PROB
    predictions = specificity_prediction_helper.predict(a_domain)
    #return {res: ((log(prob) - log(apriori_prob[res]))
    #              if prob > 0.0 else -float('inf'))
    #        for res, prob in predictions.items()}
    return {res: log(prob) if prob > 0.0 else -float('inf')
            for res, prob in predictions.items()}


def build_bgc_assembly_line(fragmented_bgc: Fragmented_BGC_Cluster,
                            specificity_prediction_helper: SpecificityPredictionHelper)\
        -> List[BGC_Module]:

    genes = list(chain(*fragmented_bgc.fragments))
    iterative_genes = get_iterative_genes(genes)
    iterative_modules_ids = get_iterative_modules_ids(genes)
    module_mods = get_modules_modifications(genes)
    module_id_to_genomic_context = get_modules_genomic_context(genes)

    modules = []
    for bgc_fragment_idx, bgc_fragment_genes in enumerate(fragmented_bgc.fragments):
        for gene_idx, gene in enumerate(bgc_fragment_genes):

            a_domain_idx = -1  # not the same as module_idx because modules with no a_domain are skipped
            for module_idx, module in enumerate(gene.modules):
                if module.a_domain is None:
                    continue
                a_domain_idx += 1
                module_id = BGC_Module_ID(gene.gene_id, module_idx)
                module_genomic_context = module_id_to_genomic_context[module_id]

                residue_scores = get_residue_scores(module.a_domain, specificity_prediction_helper)
                modules.append(BGC_Module(gene_id=gene.gene_id,
                                          fragment_idx=bgc_fragment_idx,
                                          a_domain_idx=a_domain_idx,
                                          genomic_context=module_genomic_context,
                                          modifications=module_mods[module_id],
                                          iterative_module=module_id in iterative_modules_ids,
                                          iterative_gene=gene.gene_id in iterative_genes,
                                          aa10_code=module.a_domain.aa10,
                                          aa34_code=module.a_domain.aa34,
                                          residue_score=residue_scores))

    return modules


def build_bgc_variants(bgc: BGC_Cluster,
                       specificity_prediction_helper: SpecificityPredictionHelper,
                       antismash_cfg: antiSMASH_Processing_Config,
                       log: NerpaLogger) -> List[BGC_Variant]:
    if not any(module.a_domain is not None
               for gene in bgc.genes
               for module in gene.modules):
        log.info(f'WARNING: BGC {bgc.bgc_id} has no genes with A-domains. Skipping.')
        return []
    fragmented_bgcs = split_and_reorder(bgc, antismash_cfg, log)

    return [BGC_Variant(bgc_variant_id=BGC_Variant_ID(bgc_id=bgc.bgc_id, variant_idx=idx),
                        modules=assembly_line)
            for idx, fragmented_bgc in enumerate(fragmented_bgcs)
            if (assembly_line := build_bgc_assembly_line(fragmented_bgc,
                                                         specificity_prediction_helper))]
