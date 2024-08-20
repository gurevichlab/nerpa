from typing import List, Dict, Callable
import math
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
                     config: SpecificityPredictionConfig) -> Dict[MonomerResidue, LogProb]:
    predictions = {res: math.e ** score for res, score in _predictions.items()}
    if config.apply_step_function:
        step_function = create_step_function(config.calibration_step_function_steps)
        predictions = {res: step_function(score)
                       for res, score in predictions.items()}
    if config.normalize_scores:
        total_prob = sum(predictions.values())
        if total_prob < 1:
            predictions = {res: score + (1 - total_prob) * config.apriori_residue_prob[res]
                           for res, score in predictions.items()}
        else:
            predictions = {res: score / total_prob
                           for res, score in predictions.items()}

    if config.pseudo_counts:
        predictions = {res: score * (1 - config.pseudo_count_fraction) +
                            config.pseudo_count_fraction * config.apriori_residue_prob[res]
                       for res, score in predictions.items()}

    if config.compute_evidence:
        predictions = {res: score / config.apriori_residue_prob[res]
                       for res, score in predictions.items()}

    predictions = {res: math.log(score) for res, score in predictions.items()}
    return predictions
