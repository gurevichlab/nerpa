import math
from math import log
from typing import (
    Dict,
    List,
    NamedTuple,
    Union, Optional, Set
)

import dacite

from src.data_types import (
    Chirality,
    LogProb,
    MonomerResidue, Prob,
)
from src.antismash_parsing.genomic_context import ModuleGenomicContext, ModuleGenomicContextFeature
from src.matching.hmm_auxiliary_types import DetailedHMMEdgeType
from src.matching.alignment_step_type import MatchDetailedScore

from dataclasses import dataclass
from pathlib import Path
import yaml

from src.monomer_names_helper import NRP_Monomer, MonomerNamesHelper, MonomersDefaultFrequencies

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
    monomer_detailed_default_score: Dict[NRP_Monomer, MatchDetailedScore]
    #norine_monomers_info: NorineStats


def load_methylation_score(cfg: dict) -> Dict[MethylationMatch, LogProb]:
    meth_freqs = cfg['emission_parameters']['modification_frequencies']['METHYLATION']
    return {MethylationMatch(bgc_meth=bgc_meth, nrp_meth=nrp_meth):
                log(meth_freqs[f'BGC_{bgc_meth}_NRP_{nrp_meth}'])
            for bgc_meth in (False, True)
            for nrp_meth in (False, True)}


def load_chirality_score(cfg: dict,
                         default_freqs: MonomersDefaultFrequencies) -> Dict[ChiralityMatch, LogProb]:
    epim_freqs = cfg['emission_parameters']['modification_frequencies']['EPIMERIZATION']
    default_d_chr_freq = default_freqs.d_chirality

    chr_score = {ChiralityMatch(bgc_epim=bgc_epim, nrp_chr=nrp_chr):
                     log(epim_freqs[f'BGC_{bgc_epim}_NRP_{npr_epim}'])
                 for bgc_epim in (False, True)
                 for nrp_chr, npr_epim in ((Chirality.L, False), (Chirality.D, True))}
    chr_score.update({
            ChiralityMatch(bgc_epim=bgc_epim, nrp_chr=Chirality.UNKNOWN):
                (default_d_chr_freq * log(epim_freqs[f'BGC_{bgc_epim}_NRP_True']) +
                 (1 - default_d_chr_freq) * log(epim_freqs[f'BGC_{bgc_epim}_NRP_False']))
            for bgc_epim in (False, True)
    })
    return chr_score


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


def compute_monomers_default_detailed_score(default_freqs: MonomersDefaultFrequencies) \
        -> Dict[NRP_Monomer, MatchDetailedScore]:
    monomer_default_score = {}
    for residue, residue_freq in default_freqs.residue.items():
        for methylation in (False, True):
            for chirality in Chirality:
                residue_score = log(residue_freq)
                methylation_score = log(default_freqs.methylation) \
                    if methylation else log(1 - default_freqs.methylation)

                match chirality:
                    case Chirality.D:
                        chirality_score = log(default_freqs.d_chirality)
                    case Chirality.L:
                        # slightly inaccurate because some monomers are stereosymmetric
                        chirality_score = log(1 - default_freqs.d_chirality)
                    case Chirality.UNKNOWN:
                        chirality_score = default_freqs.d_chirality * log(default_freqs.d_chirality) + \
                                          (1 - default_freqs.d_chirality) * log(1 - default_freqs.d_chirality)

                monomer = NRP_Monomer(residue=residue,
                                      methylated=methylation,
                                      chirality=chirality)
                monomer_default_score[monomer] = MatchDetailedScore(residue_score,
                                                                    methylation_score,
                                                                    chirality_score)

    # assert that all scores sum to 1
    total_score = sum(math.e ** (score.residue_score + score.methylation_score + score.chirality_score)
                      for mon, score in monomer_default_score.items()
                      if mon.chirality != Chirality.UNKNOWN)
    assert math.isclose(total_score, 1.0), \
        f'Total score of monomers default score is {total_score}, expected 1.0'

    return monomer_default_score

def load_hmm_scoring_config(path_to_config: Path,
                            monomer_names_helper: MonomerNamesHelper) -> HMMScoringConfig:
    cfg = yaml.safe_load(path_to_config.open('r'))

    default_scores = compute_monomers_default_detailed_score(monomer_names_helper.default_frequencies)
    return HMMScoringConfig(methylation_score=load_methylation_score(cfg),
                            chirality_score=load_chirality_score(cfg, monomer_names_helper.default_frequencies),
                            edge_weight_parameters=load_edge_weight_params(cfg),
                            relevant_genomic_features=EDGE_TYPE_DEPENDENCIES,
                            monomer_detailed_default_score=default_scores)
