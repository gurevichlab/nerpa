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
    SVM_Prediction, BGC_Module_ID
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
from src.antismash_parsing.genomic_context import (
    FragmentGenomicContext,
    GeneGenomicContext,
    get_bgc_fragment_genomic_context,
    get_gene_genomic_context,
    get_module_genomic_context,
)
from itertools import chain


def get_residue_scores(a_domain: A_Domain,
                       specificity_prediction_helper: SpecificityPredictionHelper) -> Dict[MonomerResidue, LogProb]:
    apriori_prob = specificity_prediction_helper.config.APRIORI_RESIDUE_PROB
    predictions = specificity_prediction_helper.predict(a_domain)
    return {res: ((log(prob) - log(apriori_prob[res]))
                  if prob > 0.0 else -float('inf'))
            for res, prob in predictions.items()}


def build_bgc_assembly_line(bgc_fragments: List[BGC_Cluster],
                            specificity_prediction_helper: SpecificityPredictionHelper)\
        -> List[BGC_Module]:

    iterative_genes = get_iterative_genes(bgc_fragments)
    iterative_modules_ids = get_iterative_modules_ids(bgc_fragments)
    module_mods = get_modules_modifications(bgc_fragments)

    modules = []
    for bgc_fragment_idx, bgc_fragment in enumerate(bgc_fragments):
        fragment_context = get_bgc_fragment_genomic_context(bgc_fragment_idx, bgc_fragments)
        for gene_idx, gene in enumerate(bgc_fragment.genes):
            gene_context = get_gene_genomic_context(gene_idx,
                                                    bgc_fragment.genes,
                                                    fragment_context)

            a_domain_idx = -1  # not the same as module_idx because modules with no a_domain are skipped
            for module_idx, module in enumerate(gene.modules):
                if module.a_domain is None:
                    continue
                a_domain_idx += 1
                module_id = BGC_Module_ID(gene.gene_id, module_idx)
                module_genomic_context = get_module_genomic_context(module_idx, gene, gene_context)

                residue_scores = get_residue_scores(module.a_domain, specificity_prediction_helper)
                modules.append(BGC_Module(gene_id=gene.gene_id,
                                          fragment_idx=bgc_fragment_idx,
                                          a_domain_idx=a_domain_idx,
                                          genomic_context=module_genomic_context,
                                          modifications=module_mods[module_id],
                                          iterative_module=a_domain_idx in iterative_modules_ids,
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
        log.info(f'WARNING: BGC {bgc.bgc_idx} has no genes with A-domains. Skipping.')
        return []
    raw_fragmented_bgcs = split_and_reorder(bgc, antismash_cfg, log)

    bgc_id = BGC_ID(genome_id=bgc.genome_id,
                    contig_idx=bgc.contig_idx,
                    bgc_idx=bgc.bgc_idx)

    return [BGC_Variant(bgc_variant_id=BGC_Variant_ID(bgc_id=bgc_id, variant_idx=idx),
                        modules=assembly_line)
            for idx, raw_fragmented_bgc in enumerate(raw_fragmented_bgcs)
            if (assembly_line := build_bgc_assembly_line(raw_fragmented_bgc,
                                                         specificity_prediction_helper))]
