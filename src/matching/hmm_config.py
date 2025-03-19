from typing import (
    Dict,
    List,
    NamedTuple,
    Union
)

from src.data_types import (
    Chirality,
    LogProb,
    MonomerResidue,
)
from src.antismash_parsing.location_features import (
    BGC_Fragment_Loc_Feature,
    GeneLocFeature,
    ModuleLocFeature,
)
from src.matching.hmm_auxiliary_types import DetailedHMMEdgeType

from dataclasses import dataclass
from pathlib import Path
import yaml


class ChiralityMatch(NamedTuple):
    bgc_epim: bool
    nrp_chr: Chirality


class MethylationMatch(NamedTuple):
    bgc_meth: bool
    nrp_meth: bool


SingleFeatureContext = Union[ModuleLocFeature, GeneLocFeature, BGC_Fragment_Loc_Feature, None]
EdgeWeightsParams = Dict[DetailedHMMEdgeType, Dict[SingleFeatureContext, float]]


@dataclass
class HMMScoringConfig:
    methylation_score: Dict[MethylationMatch, LogProb]
    chirality_score: Dict[ChiralityMatch, LogProb]
    edge_weight_parameters: EdgeWeightsParams


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


def load_edge_weight_params(cfg: dict) -> EdgeWeightsParams:
    parsed_data = {}
    for edge_type, context_to_probs in cfg['edge_weight_parameters'].items():
        parsed_data[DetailedHMMEdgeType[edge_type]] = {}
        for context_str, prob in context_to_probs.items():
            if context_str == 'None':
                parsed_data[DetailedHMMEdgeType[edge_type]][None] = prob
                continue
            if context_str in ModuleLocFeature.__members__:
                context = ModuleLocFeature[context_str]
            elif context_str in GeneLocFeature.__members__:
                context = GeneLocFeature[context_str]
            elif context_str in BGC_Fragment_Loc_Feature.__members__:
                context = BGC_Fragment_Loc_Feature[context_str]
            else:
                raise ValueError(f'Unknown context: {context_str}')
            parsed_data[DetailedHMMEdgeType[edge_type]][context] = prob

    return parsed_data


def load_hmm_scoring_config(path_to_config: Path) -> HMMScoringConfig:
    cfg = yaml.safe_load(path_to_config.open('r'))

    return HMMScoringConfig(methylation_score=load_methylation_score(cfg),
                            chirality_score=load_chirality_score(cfg),
                            edge_weight_parameters=load_edge_weight_params(cfg))
