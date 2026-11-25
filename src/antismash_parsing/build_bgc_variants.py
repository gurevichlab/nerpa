from math import log, exp
import math
from typing import (
    Dict,
    List
)
from src.aa_specificity_prediction_model.specificity_prediction_helper import SpecificityPredictionHelper
from src.antismash_parsing.antismash_parser_types import (
    A_Domain,
    BGC_Cluster,
    BGC_Module_ID, Fragmented_BGC_Cluster
)
from src.antismash_parsing.determine_modifications import (
    get_modules_modifications,
    get_iterative_genes,
    get_iterative_modules_ids
)
from src.antismash_parsing.bgc_variant_types import (
    BGC_Variant_ID,
    BGC_Variant,
    BGC_Module,
)
from src.general_type_aliases import LogProb, LogOdds, Prob
from src.monomer_names_helper.monomer_names_helper import NerpaResidue
from src.pipeline.logging.buffered_logger import BufferedLogger
from src.config import antiSMASH_Processing_Config
from src.antismash_parsing.bgcs_split_and_reorder import split_and_reorder
from src.antismash_parsing.genomic_context import get_modules_genomic_context
from itertools import chain


def get_residue_scores(
        a_domain: A_Domain,
        specificity_prediction_helper: SpecificityPredictionHelper
) -> Dict[NerpaResidue, LogProb]:
    #apriori_prob = specificity_prediction_helper.config.APRIORI_RESIDUE_PROB
    predictions = specificity_prediction_helper.predict(a_domain)
    #return {res: ((log(prob) - log(apriori_prob[res]))
    #              if prob > 0.0 else -float('inf'))
    #        for res, prob in predictions.items()}
    return {res: log(prob) if prob > 0.0 else -float('inf')
            for res, prob in predictions.items()}

def compute_prediction_uncertainty(
        modules: List[BGC_Module],
        specificity_prediction_helper: SpecificityPredictionHelper,
) -> Prob:
    r"""
    Computes the \textbf{mean normalized entropy across A-domain predictions}:

    $$
    u = \frac{1}{m}\sum_i
    \frac{H(p_i)-H_{\min}}{\log(n)-H_{\min}}
    $$
    
    where $m$ is the number of modules, $n$ is the number of residues, $H(p_i)$ is the entropy of the predicted residue distribution for module $i$, and $H_{\min}$ is the minimum possible entropy for an A domain (it's not 0 because zero probabilities are not allowed).

    The normalized entropy is in the range $[0, 1]$, where 0 indicates a very confident prediction (low uncertainty) and 1 indicates a very uncertain prediction (high uncertainty).
    """
    minimum_entropy: float = (
        specificity_prediction_helper.minimum_possible_prediction_entropy()
    )
    number_of_residues: int = len(
        specificity_prediction_helper.config.APRIORI_RESIDUE_PROB
    )
    maximum_entropy: float = log(number_of_residues)
    entropy_range: float = maximum_entropy - minimum_entropy

    assert entropy_range >= 0, (
        f"Maximum entropy ({maximum_entropy}) is less than minimum entropy ({minimum_entropy})."
    )

    normalized_entropies: List[float] = []
    for module in modules:
        probabilities: List[Prob] = [
            exp(log_probability)
            for log_probability in module.residue_score.values()
            if log_probability != -float('inf')
        ]
        entropy: float = -sum(
            probability * log(probability)
            for probability in probabilities
            if probability > 0
        )
        normalized_entropy: float = (
            entropy - minimum_entropy
        ) / entropy_range
        assert 0 <= normalized_entropy <= 1, (
            f"Normalized entropy ({normalized_entropy}) is not in the range [0, 1]."
        )
        normalized_entropies.append(normalized_entropy)
 

    return sum(normalized_entropies) / len(normalized_entropies)


def build_bgc_assembly_line(fragmented_bgc: Fragmented_BGC_Cluster,
                            specificity_prediction_helper: SpecificityPredictionHelper)\
        -> List[BGC_Module]:

    genes = list(chain(*fragmented_bgc.fragments))
    iterative_genes = get_iterative_genes(genes)
    iterative_modules_ids = get_iterative_modules_ids(genes)
    module_mods = get_modules_modifications(genes)
    module_id_to_genomic_context = get_modules_genomic_context(fragmented_bgc.fragments)

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
                       log: BufferedLogger,
                       let_it_crash: bool = False) -> List[BGC_Variant]:
    if not any(module.a_domain is not None
               for gene in bgc.genes
               for module in gene.modules):
        log.warning(f'BGC {bgc.bgc_id} has no genes with A-domains. Skipping.')
        return []

    try:
        fragmented_bgcs = split_and_reorder(bgc, antismash_cfg, log)

        bgc_variants: list[BGC_Variant] = [
            BGC_Variant(
                bgc_variant_id=BGC_Variant_ID(bgc_id=bgc.bgc_id, variant_idx=idx),
                modules=assembly_line,
                metadata=bgc.metadata,
                prediction_uncertainty=(
                    compute_prediction_uncertainty(assembly_line, specificity_prediction_helper)
                )
            )
            for idx, fragmented_bgc in enumerate(fragmented_bgcs)
            if (assembly_line := build_bgc_assembly_line(fragmented_bgc,
                                                         specificity_prediction_helper))
        ]
    except Exception as e:
        log.error(f'Error processing BGC {bgc.bgc_id}. {e.__class__.__name__}: {e}. Skipping.')
        if let_it_crash:
            raise e
        return []

    return bgc_variants
