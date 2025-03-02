from typing import (
    Any,
    Callable,
    Dict,
    List,
    NamedTuple
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
    get_modules_modifications
)
from src.data_types import (
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
from src.antismash_parsing.location_features import (
    BGC_Fragment_Loc_Features,
    GeneLocFeatures,
    get_bgc_fragment_loc_features,
    get_gene_loc_features,
    get_module_loc_features,
)
from itertools import chain
import pandas as pd


# TODO: figure out how to write it as "lambda args: args['a_domain'].aa34"
@cached_by_key(key=lambda *args, **kwargs: args[0].aa10)
def get_residue_scores(a_domain: A_Domain,
                       residue_scoring_model: Callable[[pd.DataFrame], Dict[MonomerResidue, LogProb]],
                       monomer_names_helper: MonomerNamesHelper,
                       config: SpecificityPredictionConfig) -> Dict[MonomerResidue, LogProb]:
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

    condensed_aa_codes = {residue: ([], []) for residue in monomer_names_helper.supported_residues}
    for aa_name, aa10_codes in config.KNOWN_AA10_CODES.items():
        if aa_name == 'oxoDec':
            pass
        aa34_codes = config.KNOWN_AA34_CODES[aa_name]
        residue = monomer_names_helper.parsed_name(aa_name, 'antismash').residue
        condensed_aa_codes[residue][0].extend(aa10_codes)
        condensed_aa_codes[residue][1].extend(aa34_codes)

    for residue, (aa10_codes, aa34_codes) in condensed_aa_codes.items():
        try:
            aa_code_scores = [max(similarity_score(aa_code, known_aa_code)
                                  for known_aa_code in known_aa_codes)
                              for aa_code, known_aa_codes in [(a_domain.aa10, aa10_codes),
                                                              (a_domain.aa34, aa34_codes)]]
        except:
            pass
        svm_scores = [svm_score(residue, a_domain.svm[level])
                      for level in SVM_LEVEL]
        scoring_table.loc[residue] = aa_code_scores + svm_scores

    result = residue_scoring_model(scoring_table)
    return result


def build_gene_assembly_line(gene: Gene,
                             gene_loc_features: GeneLocFeatures,
                             fragment_idx: int,
                             residue_scoring_model: Any,
                             monomer_names_helper: MonomerNamesHelper,
                             config: SpecificityPredictionConfig) -> List[BGC_Module]:
    iterative_modules_idxs = get_iterative_modules_idxs(gene)
    modules_modifications = get_modules_modifications(gene)

    built_modules = []
    a_domain_idx = -1  # indexes start from 0
    for module_idx, module in enumerate(gene.modules):
        if module.a_domain is None:
            continue
        a_domain_idx += 1

        residue_scores = get_residue_scores(module.a_domain,
                                            residue_scoring_model,
                                            monomer_names_helper,
                                            config)

        built_modules.append(BGC_Module(gene_id=gene.gene_id,
                                        fragment_idx=fragment_idx,
                                        module_loc=get_module_loc_features(module_idx, gene, gene_loc_features),
                                        a_domain_idx=a_domain_idx,
                                        residue_score=residue_scores,
                                        modifications=modules_modifications[module_idx],
                                        iterative_module=module_idx in iterative_modules_idxs,
                                        iterative_gene=gene.is_iterative,
                                        aa10_code=module.a_domain.aa10,
                                        aa34_code=module.a_domain.aa34))

    return built_modules


def build_bgc_fragment_assembly_line(bgc_genes: List[Gene],
                                     fragment_features: BGC_Fragment_Loc_Features,
                                     fragment_idx: int,
                                     residue_scoring_model: Any,
                                     monomer_names_helper: MonomerNamesHelper,
                                     config: SpecificityPredictionConfig) -> List[BGC_Module]:

    return list(chain.from_iterable(build_gene_assembly_line(gene,
                                                             get_gene_loc_features(gene_idx, bgc_genes, fragment_features),
                                                             fragment_idx,
                                                             residue_scoring_model,
                                                             monomer_names_helper,
                                                             config)
                                    for gene_idx, gene in enumerate(bgc_genes)))


def build_bgc_assembly_line(raw_fragmented_bgc: List[BGC_Cluster],
                            residue_scoring_model: Any,
                            monomer_names_helper: MonomerNamesHelper,
                            config: SpecificityPredictionConfig) -> List[BGC_Module]:

    return list(chain(*(build_bgc_fragment_assembly_line(bgc_fragment.genes,
                                                         fragment_features=get_bgc_fragment_loc_features(fgmnt_idx, raw_fragmented_bgc),
                                                         fragment_idx=fgmnt_idx,
                                                         residue_scoring_model=residue_scoring_model,
                                                         monomer_names_helper=monomer_names_helper,
                                                         config=config)
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
                       residue_scoring_model: Any,
                       monomer_names_helper: MonomerNamesHelper,
                       as_cfg: antiSMASH_Processing_Config,
                       spec_pred_cfg: SpecificityPredictionConfig,
                       log: NerpaLogger) -> List[BGC_Variant]:  # TODO: replace Any with proper type
    if not any(module.a_domain is not None
               for gene in bgc.genes
               for module in gene.modules):
        log.info(f'WARNING: BGC {bgc.bgc_idx} has no genes with A-domains. Skipping.')
        return []
    raw_fragmented_bgcs = split_and_reorder(bgc, as_cfg, log)
    return [BGC_Variant(genome_id=bgc.genome_id,
                        contig_idx=bgc.contig_idx,
                        bgc_idx=bgc.bgc_idx,
                        variant_idx=idx,
                        modules=assembly_line)
            for idx, raw_fragmented_bgc in enumerate(raw_fragmented_bgcs)
            if (assembly_line := build_bgc_assembly_line(raw_fragmented_bgc,
                                                         residue_scoring_model,
                                                         monomer_names_helper,
                                                         spec_pred_cfg))]
