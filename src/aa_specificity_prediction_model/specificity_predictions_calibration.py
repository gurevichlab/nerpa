from typing import (
    Callable,
    Dict,
    List,
    Literal
)
from math import log, e
from src.config import SpecificityPredictionConfig
from src.data_types import LogProb
from src.monomer_names_helper import MonomerResidue


def create_step_function(steps: List[float]) -> Callable[[float], float]:
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


def calibrate_scores(_predictions: Dict[MonomerResidue, LogProb],
                     config: SpecificityPredictionConfig,
                     model: Literal['paras', 'nerpa']) -> Dict[MonomerResidue, LogProb]:
    # convert log-probabilities to probabilities
    predictions = {res: e ** score for res, score in _predictions.items()}

    step_function_steps = config.CALIBRATION_STEP_FUNCTION_STEPS_PARAS \
        if model == 'paras' else config.CALIBRATION_STEP_FUNCTION_STEPS_NERPA
    step_function = create_step_function(step_function_steps)

    # calibrate scores to better represent the probability of residue incorporation
    predictions = {res: step_function(score)
                   for res, score in predictions.items()}

    # normalize probabilities
    total_prob = sum(predictions.values())
    if total_prob < 1:
        predictions = {res: score + (1 - total_prob) * config.APRIORI_RESIDUE_PROB[res]
                       for res, score in predictions.items()}
    else:
        predictions = {res: score / total_prob
                       for res, score in predictions.items()}

    # add pseudo-counts to avoid (near) zero probabilities
    predictions = {res: score * (1 - config.PSEUDO_COUNT_FRACTION) +
                        config.PSEUDO_COUNT_FRACTION * config.APRIORI_RESIDUE_PROB[res]
                   for res, score in predictions.items()}

    # compute log-evidence --- how likely is residue to be incorporated by the A domain compared to random chance
    predictions = {res: log(score) - log(config.APRIORI_RESIDUE_PROB[res])
                   for res, score in predictions.items()}

    return predictions
