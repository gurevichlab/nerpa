from collections import defaultdict
from pathlib import Path
from typing import Literal, Dict, List, Callable, Optional

import pandas as pd

from src.aa_specificity_prediction_model.model_wrapper import ModelWrapper
from src.antismash_parsing.antismash_name_mappings import KNOWN_SUBSTRATES
from src.antismash_parsing.antismash_parser_types import A_Domain, SVM_Prediction, SVM_LEVEL
from src.config import SpecificityPredictionConfig
from src.data_types import LogProb, Prob
from src.generic.string import hamming_distance
from src.monomer_names_helper import AA34, MonomerResidue, MonomerNamesHelper, paras_residue_to_nerpa_residue
from src.paras.paras_wrapper import ParasWrapper
from src.pipeline.logger import NerpaLogger
from src.pipeline.paras_parsing import PARAS_RESIDUE


def create_step_function(steps: List[float]) -> Callable[[Prob], float]:
    '''
     Create a step function [0,1] -> R given the list of steps
    '''
    step_len = 1 / len(steps)

    def step_function(x: float) -> float:
        for i, step in enumerate(steps):
            if x < (i + 1) * step_len:
                return steps[i]
        return steps[-1]
    return step_function


class SpecificityPredictionHelper:
    config: SpecificityPredictionConfig
    KNOWN_SPECIFICITIES: Dict[AA34, List[MonomerResidue]]
    external_predictions: Dict[AA34, Dict[MonomerResidue, LogProb]]
    DEFAULT_MODEL: Literal['nerpa', 'paras']
    log: NerpaLogger

    _calibration_step_function_nerpa: Callable[[float], float]
    _calibration_step_function_paras: Callable[[float], float]

    _cache: Dict[AA34, Dict[MonomerResidue, LogProb]]
    nerpa_model_wrapper: Optional[ModelWrapper] = None
    paras_model_wrapper: Optional[ParasWrapper] = None


    def __init__(self,
                 config: SpecificityPredictionConfig,
                 monomer_names_helper: MonomerNamesHelper,
                 log: NerpaLogger,
                 external_predictions: Optional[Dict[AA34, Dict[MonomerResidue, Prob]]] = None):
        self.config = config
        self.monomer_names_helper = monomer_names_helper
        self.external_predictions = external_predictions\
            if external_predictions is not None else {}

        self.KNOWN_SPECIFICITIES = defaultdict(list)
        for res, aa34_codes in config.KNOWN_AA34_CODES.items():
            for aa34_code in aa34_codes:
                self.KNOWN_SPECIFICITIES[aa34_code].append(res)

        if config.paras_model is not None:
            self.DEFAULT_MODEL = 'paras'
            self.paras_model_wrapper = ParasWrapper(config.paras_model)
        else:
            self.DEFAULT_MODEL = 'nerpa'
            self.nerpa_model_wrapper = ModelWrapper(config.nerpa_specificity_prediction_model)

        self._calibration_step_function_paras = create_step_function(config.CALIBRATION_STEP_FUNCTION_STEPS_PARAS)
        self._calibration_step_function_nerpa = create_step_function(config.CALIBRATION_STEP_FUNCTION_STEPS_NERPA)

        self._cache = {}
        self.log = log

    def predict(self,
                a_domain: A_Domain,
                _no_cache: bool = False,  # for debugging/training purposes
                _no_calibration: bool = False) -> Dict[MonomerResidue, Prob]:
        if a_domain.aa34 in self._cache and not _no_cache:
            return self._cache[a_domain.aa34]

        if self.config.ENABLE_DICTIONARY_LOOKUP and a_domain.aa34 in self.KNOWN_SPECIFICITIES:
            raw_predictions = {res: (1.0 if res in self.KNOWN_SPECIFICITIES[a_domain.aa34] else 0.0)
                               for res in self.monomer_names_helper.supported_residues}
            calibration_function = self._calibration_step_function_paras \
                if self.DEFAULT_MODEL == 'paras' else self._calibration_step_function_nerpa
        elif a_domain.aa34 in self.external_predictions:
            raw_predictions = self.external_predictions[a_domain.aa34]
            calibration_function = self._calibration_step_function_paras
        else:
            if self.external_predictions:
                self.log.info(f"Warning: no external prediction for {a_domain.aa34}"
                         f" - using {self.DEFAULT_MODEL} model instead.")
            if self.DEFAULT_MODEL == 'paras':
                raw_predictions = self._predict_paras(a_domain)
                calibration_function = self._calibration_step_function_paras
            else:
                raw_predictions = self._predict_nerpa(a_domain)
                calibration_function = self._calibration_step_function_nerpa

        if self.config.ENABLE_CALIBRATION and not _no_calibration:
            predictions = self.calibrate_scores(raw_predictions, calibration_function)
        else:
            predictions = raw_predictions

        if not _no_cache:
            self._cache[a_domain.aa34] = predictions
        return predictions

    def _predict_nerpa(self, a_domain: A_Domain) -> Dict[MonomerResidue, Prob]:
        def svm_score(aa_name: MonomerResidue, svm_level_prediction: SVM_Prediction) -> float:
            svm_prediction_substrates = {self.monomer_names_helper.parsed_name(monomer_name, 'antismash').residue
                                         for monomer_name in svm_level_prediction.substrates}
            if aa_name not in self.config.SVM_SUBSTRATES:
                return self.config.SVM_NOT_SUPPORTED_SCORE
            if aa_name in svm_prediction_substrates:
                return svm_level_prediction.score
            else:
                return self.config.SVM_NO_PREDICTION_SCORE

        def similarity_score(aa_code1: str, aa_code2: str) -> float:
            return 1.0 - hamming_distance(aa_code1, aa_code2) / len(aa_code1)

        scoring_table = (pd.DataFrame([], columns=self.config.SCORING_TABLE_COLUMNS)
                         .set_index(self.config.SCORING_TABLE_INDEX))

        for aa_name, aa10_codes in self.config.KNOWN_AA10_CODES.items():
            aa34_codes = self.config.KNOWN_AA34_CODES[aa_name]

            aa_code_scores = [max(similarity_score(aa_code, known_aa_code)
                                  for known_aa_code in known_aa_codes)  # type: ignore
                              for aa_code, known_aa_codes in [(a_domain.aa10, aa10_codes),
                                                              (a_domain.aa34, aa34_codes)]]
            svm_scores = [svm_score(aa_name, a_domain.svm[level])
                          for level in SVM_LEVEL]
            scoring_table.loc[aa_name] = aa_code_scores + svm_scores

        predictions = self.nerpa_model_wrapper(scoring_table, self.monomer_names_helper)
        return predictions

    def _predict_paras(self, a_domain: A_Domain) -> Dict[MonomerResidue, Prob]:
        paras_predictions = self.paras_model_wrapper.predict(a_domain.aa34)
        return self.paras_predictions_to_nerpa_predictions(paras_predictions)


    def calibrate_scores(self,
                         _predictions: Dict[MonomerResidue, Prob],
                         calibration_function: Callable[[Prob], Prob])\
            -> Dict[MonomerResidue, Prob]:
        # calibrate scores to better represent the probability of residue incorporation
        predictions = {res: calibration_function(score)
                       for res, score in _predictions.items()}

        # normalize probabilities
        total_prob = sum(predictions.values())
        if total_prob < 1:
            predictions = {res: score + (1 - total_prob) * self.config.APRIORI_RESIDUE_PROB[res]
                           for res, score in predictions.items()}
        else:
            predictions = {res: score / total_prob
                           for res, score in predictions.items()}

        # add pseudo-counts to avoid (near) zero probabilities
        predictions = {res: score * (1 - self.config.PSEUDO_COUNT_FRACTION) +
                            self.config.PSEUDO_COUNT_FRACTION * self.config.APRIORI_RESIDUE_PROB[res]
                       for res, score in predictions.items()}

        return predictions

    def paras_predictions_to_nerpa_predictions(self,
                                               paras_predictions: Dict[PARAS_RESIDUE, Prob]) -> Dict[MonomerResidue, Prob]:
        nerpa_predictions = {res: 0.0 for res in self.monomer_names_helper.supported_residues}

        for paras_residue, prob in paras_predictions.items():
            nerpa_res = paras_residue_to_nerpa_residue(paras_residue, self.monomer_names_helper)
            nerpa_predictions[nerpa_res] += prob

        for res in nerpa_predictions:
            if nerpa_predictions[res] > 1:
                print(f"Invalid probability: {nerpa_predictions[res]}")
                nerpa_predictions[res] = 1.0

        return nerpa_predictions
