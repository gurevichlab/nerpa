import pandas as pd
import numpy as np
from pathlib import Path
import joblib

from src.data_types import ResidueScores


class ModelWrapper(object):
    def __init__(self, 
                 model_dump: Path, 
                 lookup_score: str = "aa34_score",
                 lookup_threshold: float = 1e-3):
        """Wrapper over sklearn.ensemble.RandomForestClassifier that combines lookup with the prediction of log_probs

        Parameters
        ----------
        model_dump : Path
            Path to joblib.dump'ed RandomForestClassifier model.
        lookup_score : {"aa34_score", "aa10_score"}, optional
            Score to use for the lookup. Defaults to "aa34_score".
        lookup_threshold : float, optional
            Allowed error for the lookup_score. Residues with `score > (1-thres)` will be assumed to be correct matches
            and assigned `prob = 1.0`. Defaults to 0.001 which requires 100% match for both scores.
        """
        self.model = joblib.load(model_dump)
        self.lookup_col = lookup_score
        self.lookup_threshold = lookup_threshold
        self.calibration_probs = np.array([0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
            0.001, 0.001, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 
            0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 
            0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 
            0.003, 0.003, 0.005, 0.005, 0.005, 0.005, 0.005, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.097, 0.097, 0.097, 
            0.097, 0.097, 0.097, 0.097, 0.097, 0.097, 0.104, 0.104, 0.167, 0.167, 0.167, 0.167, 0.167, 0.167, 0.167, 0.167, 
            0.167, 0.167, 0.167, 0.167, 0.877, 0.877, 0.877])
        self.calibration_bins = np.linspace(0, 1, self.calibration_probs.shape[0])

    def __call__(self, scoring_table: pd.DataFrame) -> ResidueScores:
        scores = pd.Series(0, index=scoring_table.index, dtype=float)
        
        # Filter rows with no match in the lookup table
        mask_predict = (1 - scoring_table[self.lookup_col]).abs() > self.lookup_threshold

        # Get predictions from the model keeping only the pos class prediction
        preds = self.model.predict_proba(scoring_table.loc[mask_predict])[:,1]
        
        # Calibrate probabilities
        bin_ids = np.digitize(preds, self.calibration_bins, right=True)
        scores.loc[mask_predict] = np.log(self.calibration_probs[bin_ids])

        # Convert to ResidueScores=Dict[str, float] and return
        return scores.to_dict()

