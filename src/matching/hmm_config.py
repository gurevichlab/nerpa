from typing import (
    Dict,
    List,
    NamedTuple,
    Union, Optional, Set
)

from src.data_types import (
    Chirality,
    LogProb,
    MonomerResidue,
)
from src.antismash_parsing.genomic_context import ModuleGenomicContext, ModuleGenomicContextFeature
from src.matching.hmm_auxiliary_types import DetailedHMMEdgeType

from dataclasses import dataclass
from pathlib import Path
import yaml


# TODO: put this in config
EDGE_TYPE_DEPENDENCIES = {
    # Insertions
    DetailedHMMEdgeType.START_INSERTING_AT_START: {
        ModuleGenomicContextFeature.PKS_UPSTREAM,
    },
    DetailedHMMEdgeType.INSERT_AT_START: set(),
    DetailedHMMEdgeType.START_INSERTING: {
        ModuleGenomicContextFeature.END_OF_GENE,
        ModuleGenomicContextFeature.PKS_DOWNSTREAM,
    },
    DetailedHMMEdgeType.INSERT: set(),
    DetailedHMMEdgeType.START_INSERTING_AT_END: {ModuleGenomicContextFeature.PKS_DOWNSTREAM},
    DetailedHMMEdgeType.INSERT_AT_END: set(),

    # Skips
    DetailedHMMEdgeType.START_SKIP_MODULES_AT_START: set(),
    DetailedHMMEdgeType.SKIP_MODULE_AT_START: {ModuleGenomicContextFeature.ONLY_A_DOMAIN},
    DetailedHMMEdgeType.SKIP_MODULE: {ModuleGenomicContextFeature.ONLY_A_DOMAIN},
    DetailedHMMEdgeType.SKIP_GENE: set(),

    DetailedHMMEdgeType.END_MATCHING: set(),
    DetailedHMMEdgeType.SKIP_MODULE_END_MATCHING: set(),
    DetailedHMMEdgeType.SKIP_MODULE_AT_END: {ModuleGenomicContextFeature.ONLY_A_DOMAIN},

    # Iterations
    DetailedHMMEdgeType.ITERATE_MODULE: set(),
    DetailedHMMEdgeType.ITERATE_GENE: set(),

    # Other (not needed for training)
    #DetailedHMMEdgeType.MATCH: {ModuleGenomicContextFeature.ONLY_A_DOMAIN},
    #DetailedHMMEdgeType.START_MATCHING: set(),
    #DetailedHMMEdgeType.NO_INSERTIONS: set(),
    #DetailedHMMEdgeType.END_INSERTING: set(),
}


class ChiralityMatch(NamedTuple):
    bgc_epim: bool
    nrp_chr: Chirality


class MethylationMatch(NamedTuple):
    bgc_meth: bool
    nrp_meth: bool


@dataclass
class HMMScoringConfig:
    methylation_score: Dict[MethylationMatch, LogProb]
    chirality_score: Dict[ChiralityMatch, LogProb]
    edge_weight_parameters: Dict[DetailedHMMEdgeType, Dict[ModuleGenomicContext, LogProb]]
    relevant_genomic_features: Dict[DetailedHMMEdgeType, Set[ModuleGenomicContextFeature]]


def load_methylation_score(cfg: dict) -> Dict[MethylationMatch, LogProb]:
    return {MethylationMatch(bgc_meth=bgc_meth, nrp_meth=nrp_meth):
                cfg['emission_parameters']['methylation_score'][f'BGC_{bgc_meth}_NRP_{nrp_meth}']
            for bgc_meth in (False, True)
            for nrp_meth in (False, True)}


def load_chirality_score(cfg: dict) -> Dict[ChiralityMatch, LogProb]:
    return {ChiralityMatch(bgc_epim=bgc_epim, nrp_chr=nrp_chr):
                cfg['emission_parameters']['chirality_score'][f'BGC_{bgc_epim}_NRP_{nrp_chr.name}']
            for bgc_epim in (False, True)
            for nrp_chr in Chirality}


def load_edge_weight_params(cfg: dict) -> Dict[DetailedHMMEdgeType, Dict[ModuleGenomicContext, LogProb]]:
    edge_weights_cfg = cfg['edge_weight_parameters']
    ET = DetailedHMMEdgeType
    MGF = ModuleGenomicContextFeature

    parsed_data = {}
    for edge_type_name, weight_for_gc in edge_weights_cfg.items():
        parsed_data[ET[edge_type_name]] = {}
        for gc_str, weight in weight_for_gc.items():
            gc = tuple(MGF[feature_str] for feature_str in eval(gc_str))
            parsed_data[ET[edge_type_name]][gc] = weight

    return parsed_data


def load_hmm_scoring_config(path_to_config: Path) -> HMMScoringConfig:
    cfg = yaml.safe_load(path_to_config.open('r'))

    return HMMScoringConfig(methylation_score=load_methylation_score(cfg),
                            chirality_score=load_chirality_score(cfg),
                            edge_weight_parameters=load_edge_weight_params(cfg),
                            relevant_genomic_features=EDGE_TYPE_DEPENDENCIES)
