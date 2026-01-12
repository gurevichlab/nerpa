import math
from math import log
from typing import (
    Dict,
    NamedTuple
)

import pandas as pd

from src.aa_specificity_prediction_model.specificity_prediction_helper import SpecificityPredictionHelper
from src.general_type_aliases import (
    LogProb,
     Prob,
)
from src.antismash_parsing.genomic_context import ModuleGenomicContext, ModuleGenomicContextFeature
from src.generic.numeric import safe_log
from src.hmm.hmm_auxiliary_types import DetailedHMMEdgeType, DetailedHMMStateType
from src.matching.alignment_step_type import MatchDetailedScore

from dataclasses import dataclass
from pathlib import Path
import yaml

from src.monomer_names_helper import (
    Chirality,
    NRP_Monomer,
    NerpaResidue,
    MonomerNamesHelper,
    MonomersDefaultFrequencies
)

ST = DetailedHMMStateType
ET = DetailedHMMEdgeType
MGF = ModuleGenomicContextFeature

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

    # the score of a monomer averaged over all BGC modules from MIBiG
    # used for normalization
    monomer_detailed_default_score: Dict[NRP_Monomer, MatchDetailedScore]
    # probability that an inserted monomer has unknown residue
    insert_unknown_freq: Prob
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


def load_edge_weight_params(cfg: dict) -> Dict[ST, Dict[ModuleGenomicContext, Dict[ET, LogProb]]]:
    edge_weights_cfg = cfg['edge_weight_parameters']
    ET = DetailedHMMEdgeType
    MGF = ModuleGenomicContextFeature

    parsed_data = {}
    for state_type_name, ctxt_dict in edge_weights_cfg.items():
        parsed_data[ST[state_type_name]] = {}
        for gc_str, edge_types_dict in ctxt_dict.items():
            gc = tuple(MGF[feature_str] for feature_str in eval(gc_str))
            parsed_data[ST[state_type_name]][gc] = {}
            for edge_type_name, log_prob in edge_types_dict.items():
                parsed_data[ST[state_type_name]][gc][ET[edge_type_name]] = log_prob

    return parsed_data


def compute_monomers_default_detailed_score_old(default_freqs: MonomersDefaultFrequencies) \
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


def compute_monomers_default_detailed_score(paras_default_residue_freqs: Dict[NerpaResidue, Prob],
                                            methylation_score: Dict[MethylationMatch, LogProb],
                                            chirality_score: Dict[ChiralityMatch, LogProb],
                                            methylation_freq: Prob,
                                            epimerization_freq: Prob) \
        -> Dict[NRP_Monomer, MatchDetailedScore]:
    monomer_default_score = {}
    for residue, residue_freq in paras_default_residue_freqs.items():
        for nrp_meth in (False, True):
            for nrp_chr in Chirality:
                residue_score = safe_log(residue_freq)
                default_methylation_score = (
                    methylation_freq * methylation_score[MethylationMatch(bgc_meth=True, nrp_meth=nrp_meth)]
                    + (1 - methylation_freq) * methylation_score[MethylationMatch(bgc_meth=False, nrp_meth=nrp_meth)]
                )
                default_chirality_score = (
                    epimerization_freq * chirality_score[ChiralityMatch(bgc_epim=True, nrp_chr=nrp_chr)]
                    + (1 - epimerization_freq) * chirality_score[ChiralityMatch(bgc_epim=False, nrp_chr=nrp_chr)]
                )

                monomer = NRP_Monomer(residue=residue,
                                      methylated=nrp_meth,
                                      chirality=nrp_chr)
                monomer_default_score[monomer] = MatchDetailedScore(residue_score,
                                                                    default_methylation_score,
                                                                    default_chirality_score)

    return monomer_default_score


def load_paras_default_residue_scores(paras_results_tsv: Path,
                                      specificity_prediction_helper: SpecificityPredictionHelper,
                                      monomer_names_helper: MonomerNamesHelper) -> Dict[NerpaResidue, Prob]:
    assert specificity_prediction_helper.DEFAULT_MODEL == 'paras', \
        'PARAS residue frequencies can be loaded only if PARAS is used as the default specificity prediction model.'
    paras_results = pd.read_csv(paras_results_tsv, sep='\t')
    paras_predictions = [
        dict(zip(group["substrate"], group["score"]))
        for _, group in paras_results.groupby("id")
    ]

    predictions = [specificity_prediction_helper.paras_predictions_to_nerpa_predictions(single_prediction)
                   for single_prediction in paras_predictions]

    calibrated_predictions = [
        specificity_prediction_helper.calibrate_scores(prediction,
                                                       specificity_prediction_helper._calibration_step_function_paras)
        for prediction in predictions
    ]

    residues = calibrated_predictions[0].keys()
    averaged_predictions = {
        res: sum(pred[res] for pred in calibrated_predictions) / len(calibrated_predictions)
        for res in residues
    }

    '''
    print('Loaded PARAS default residue frequencies:\n',
          '\n'.join(f"{res}: {freq:.4f}"
                    for res, freq in averaged_predictions.items()))
    '''

    assert math.isclose(sum(averaged_predictions.values()), 1.0, rel_tol=1e-4), \
        f'Total score of averaged PARAS predictions is {sum(averaged_predictions.values())}, expected 1.0'

    return averaged_predictions


def load_hmm_scoring_config(nerpa_dir: Path,
                            path_to_config: Path,
                            specificity_prediction_helper: SpecificityPredictionHelper,
                            monomer_names_helper: MonomerNamesHelper) -> HMMScoringConfig:
    cfg = yaml.safe_load(path_to_config.open('r'))

    methylation_score = load_methylation_score(cfg)
    chirality_score = load_chirality_score(cfg, monomer_names_helper.default_frequencies)

    paras_prediction_table = nerpa_dir / cfg['paras_predictions_table']
    paras_default_residue_scores = load_paras_default_residue_scores(paras_prediction_table,
                                                                     specificity_prediction_helper,
                                                                     monomer_names_helper)
    default_scores = compute_monomers_default_detailed_score(paras_default_residue_freqs=paras_default_residue_scores,
                                                             methylation_freq=cfg['mibig_methylation_freq'],
                                                             epimerization_freq=cfg['mibig_epimerization_freq'],
                                                             methylation_score=methylation_score,
                                                             chirality_score=chirality_score)

    return HMMScoringConfig(methylation_score=methylation_score,
                            chirality_score=chirality_score,
                            edge_weight_parameters=load_edge_weight_params(cfg),
                            monomer_detailed_default_score=default_scores,
                            insert_unknown_freq=cfg['insert_unknown_freq'],)
