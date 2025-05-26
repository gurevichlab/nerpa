from math import e
from typing import Dict
import pandas as pd
from pathlib import Path
import joblib

from src.data_types import Prob
from src.monomer_names_helper import (
    MonomerNamesHelper,
    MonomerResidue
)


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

    def __call__(self,
                 scoring_table: pd.DataFrame,
                 monomer_names_helper: MonomerNamesHelper) -> Dict[MonomerResidue, Prob]:
        scores = pd.Series(0, index=scoring_table.index, dtype=float)
        
        # Filter rows with no match in the lookup table
        # WITH LOOKUP: mask_predict = (1 - scoring_table[self.lookup_col]).abs() > self.lookup_threshold

        # Get predictions from the model keeping only the pos class prediction
        # WITH LOOKUP: scores.loc[mask_predict] = self.model.predict_log_proba(scoring_table.loc[mask_predict])[:,1]
        scores[:] = self.model.predict_log_proba(scoring_table)[:,1]

        # Convert to ResidueScores=Dict[str, float]
        scores_dict = scores.to_dict()

        # Condense predictions to the set of supported residues
        # I take the max score for each residue, as predictions are independent
        scores = {res: 0.0 for res in monomer_names_helper.supported_residues}
        for monomer_name, log_prob in scores_dict.items():
            monomer_res = monomer_names_helper.parsed_name(monomer_name, name_format='antismash').residue
            scores[monomer_res] = max(scores[monomer_res], e ** log_prob)

        return scores

