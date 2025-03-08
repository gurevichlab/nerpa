from typing import (
    Any,
    Callable,
    Dict,
    List,
    NamedTuple,
    Optional
)
from src.aa_specificity_prediction_model.model_wrapper import ModelWrapper
from src.aa_specificity_prediction_model.specificity_predictions_calibration import calibrate_scores
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
    get_modules_modifications
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
    get_bgc_fragment_loc_features,
    get_gene_loc_features,
    get_module_loc_features,
)
from itertools import chain
import pandas as pd


# TODO: figure out how to write it as "lambda args: args['a_domain'].aa34"
@cached_by_key(key=lambda *args, **kwargs: args[0].aa34)
def get_residue_scores(a_domain: A_Domain,
                       external_specificity_predictions: Optional[Dict[AA34, Dict[MonomerResidue, LogProb]]],
                       residue_scoring_model: ModelWrapper,
                       monomer_names_helper: MonomerNamesHelper,
                       config: SpecificityPredictionConfig,
                       log: NerpaLogger) -> Dict[MonomerResidue, LogProb]:
    if external_specificity_predictions is not None:
        if a_domain.aa34 in external_specificity_predictions:
            return calibrate_scores(external_specificity_predictions[a_domain.aa34],
                                    config,
                                    model='paras')
        else:
            log.info(f'WARNING: No external specificity predictions for AA34 {a_domain.aa34}. '
                     f'Using Nerpa model.')

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

    predictions = residue_scoring_model(scoring_table,
                                        dictionary_lookup=config.ENABLE_DICTIONARY_LOOKUP,
                                        monomer_names_helper=monomer_names_helper)
    return calibrate_scores(predictions, config, model='nerpa')


def build_gene_assembly_line(gene: Gene,
                             gene_loc_features: GeneGenomicContext,
                             fragment_idx: int,
                             external_specificity_predictions: Optional[Dict[AA34, Dict[MonomerResidue, LogProb]]],
                             residue_scoring_model: ModelWrapper,
                             monomer_names_helper: MonomerNamesHelper,
                             config: SpecificityPredictionConfig,
                             log: NerpaLogger) -> List[BGC_Module]:
    iterative_modules_idxs = get_iterative_modules_idxs(gene)
    modules_modifications = get_modules_modifications(gene)

    built_modules = []
    a_domain_idx = -1  # indexes start from 0
    for module_idx, module in enumerate(gene.modules):
        if module.a_domain is None:
            continue
        a_domain_idx += 1

        residue_scores = get_residue_scores(module.a_domain,
                                            external_specificity_predictions,
                                            residue_scoring_model,
                                            monomer_names_helper,
                                            config,
                                            log)

        built_modules.append(BGC_Module(gene_id=gene.gene_id,
                                        fragment_idx=fragment_idx,
                                        genomic_context=get_module_loc_features(module_idx, gene, gene_loc_features),
                                        a_domain_idx=a_domain_idx,
                                        residue_score=residue_scores,
                                        modifications=modules_modifications[module_idx],
                                        iterative_module=module_idx in iterative_modules_idxs,
                                        iterative_gene=gene.is_iterative,
                                        aa10_code=module.a_domain.aa10,
                                        aa34_code=module.a_domain.aa34))

    return built_modules


def build_bgc_fragment_assembly_line(bgc_genes: List[Gene],
                                     fragment_features: FragmentGenomicContext,
                                     fragment_idx: int,
                                     external_specificity_predictions: Optional[Dict[AA34, Dict[MonomerResidue, LogProb]]],
                                     residue_scoring_model: ModelWrapper,
                                     monomer_names_helper: MonomerNamesHelper,
                                     config: SpecificityPredictionConfig,
                                     log: NerpaLogger) -> List[BGC_Module]:

    return list(chain.from_iterable(build_gene_assembly_line(gene,
                                                             get_gene_loc_features(gene_idx, bgc_genes, fragment_features),
                                                             fragment_idx,
                                                             external_specificity_predictions,
                                                             residue_scoring_model,
                                                             monomer_names_helper,
                                                             config,
                                                             log)
                                    for gene_idx, gene in enumerate(bgc_genes)))


def build_bgc_assembly_line(raw_fragmented_bgc: List[BGC_Cluster],
                            external_specificity_predictions: Optional[Dict[AA34, Dict[MonomerResidue, LogProb]]],
                            residue_scoring_model: ModelWrapper,
                            monomer_names_helper: MonomerNamesHelper,
                            config: SpecificityPredictionConfig,
                            log: NerpaLogger) -> List[BGC_Module]:

    return list(chain(*(build_bgc_fragment_assembly_line(bgc_fragment.genes,
                                                         fragment_features=get_bgc_fragment_loc_features(fgmnt_idx, raw_fragmented_bgc),
                                                         fragment_idx=fgmnt_idx,
                                                         external_specificity_predictions=external_specificity_predictions,
                                                         residue_scoring_model=residue_scoring_model,
                                                         monomer_names_helper=monomer_names_helper,
                                                         config=config,
                                                         log=log)
                        for fgmnt_idx, bgc_fragment in enumerate(raw_fragmented_bgc))))


# TODO: sometimes a module is split between genes (see BGC0002484). Maybe it's better to merge such modules?
def remove_genes_with_no_a_domains(bgc: BGC_Cluster) -> BGC_Cluster:
    def has_a_domains(gene: Gene) -> bool:
        return any(module.a_domain is not None for module in gene.modules)
    return BGC_Cluster(genome_id=bgc.genome_id,
                       contig_idx=bgc.contig_idx,
                       bgc_idx=bgc.bgc_idx,
                       genes=[gene for gene in bgc.genes if has_a_domains(gene)])


def build_bgc_variants(bgc: BGC_Cluster,
                       external_specificity_predictions: Optional[Dict[AA34, Dict[MonomerResidue, LogProb]]],
                       residue_scoring_model: ModelWrapper,
                       monomer_names_helper: MonomerNamesHelper,
                       antismash_cfg: antiSMASH_Processing_Config,
                       specificity_prediction_cfg: SpecificityPredictionConfig,
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
                                                         external_specificity_predictions,
                                                         residue_scoring_model,
                                                         monomer_names_helper,
                                                         specificity_prediction_cfg,
                                                         log))]
