# TODO: add licensing information


from typing import Dict, List, Sequence, Any
from pathlib import Path
import joblib
import numpy as np
from src.data_types import AA34, Prob
from src.pipeline.paras_parsing import PARAS_RESIDUE


class ParasFeaturizer(object):
    aa_params: Dict[str, Sequence[float]] = {
        "-": [0.00, 0.00, 0.00, 1, 8.3, 0.21, 13.59, 145.2, 1.00, 1.03, 0.99, 6.03, 0.06, 0.00, 0.10],
        "A": [0.07, -1.73, 0.09, 0, 8.1, -0.06, 0.00, 90.0, 1.42, 0.83, 0.66, 6.00, 0.06, -0.25, 0.25],
        "C": [0.71, -0.97, 4.13, 0, 5.5, 1.36, 1.48, 103.3, 0.70, 1.19, 1.19, 5.05, -0.56, -0.40, -0.14],
        "D": [3.64, 1.13, 2.36, 1, 13.0, -0.80, 49.70, 117.3, 1.01, 0.54, 1.46, 2.77, 0.97, -0.08, 0.08],
        "E": [3.08, 0.39, -0.07, 1, 12.3, -0.77, 49.90, 142.2, 1.51, 0.37, 0.74, 3.22, 0.85, -0.10, -0.05],
        "F": [-4.92, 1.30, 0.45, 0, 5.2, 1.27, 0.35, 191.9, 1.13, 1.38, 0.60, 5.48, -0.99, 0.18, 0.15],
        "G": [2.23, -5.36, 0.30, 0, 9.0, -0.41, 0.00, 64.9, 0.57, 0.75, 1.56, 5.97, 0.32, -0.32, 0.28],
        "H": [2.41, 1.74, 1.11, 1, 10.4, 0.49, 51.60, 160.0, 1.00, 0.87, 0.95, 7.59, 0.15, -0.03, -0.10],
        "I": [-4.44, -1.68, -1.03, 0, 5.2, 1.31, 0.13, 163.9, 1.08, 1.60, 0.47, 6.02, -1.00, -0.03, 0.10],
        "K": [2.84, 1.41, -3.14, 2, 11.3, -1.18, 49.50, 167.3, 1.16, 0.74, 1.01, 9.74, 1.00, 0.32, 0.11],
        "L": [-4.19, -1.03, -0.98, 0, 4.9, 1.21, 0.13, 164.0, 1.21, 1.30, 0.59, 5.98, -0.83, 0.05, 0.01],
        "M": [-2.49, -0.27, -0.41, 0, 5.7, 1.27, 1.43, 167.0, 1.45, 1.05, 0.60, 5.74, -0.68, -0.01, 0.04],
        "N": [3.22, 1.45, 0.84, 2, 11.6, -0.48, 3.38, 124.7, 0.67, 0.89, 1.56, 5.41, 0.70, -0.06, 0.17],
        "P": [-1.22, 0.88, 2.23, 0, 8.0, 1.1, 1.58, 122.9, 0.57, 0.55, 1.52, 6.30, 0.45, 0.23, 0.41],
        "Q": [2.18, 0.53, -1.14, 2, 10.5, -0.73, 3.53, 149.4, 1.11, 1.10, 0.98, 5.65, 0.71, -0.02, 0.12],
        "R": [2.88, 2.52, -3.44, 4, 10.5, -0.84, 52.00, 194.0, 0.98, 0.93, 0.95, 10.76, 0.80, 0.19, -0.41],
        "S": [1.96, -1.63, 0.57, 1, 9.2, -0.50, 1.67, 95.4, 0.77, 0.75, 1.43, 5.68, 0.48, -0.15, 0.23],
        "T": [0.92, -2.09, -1.40, 1, 8.6, -0.27, 1.66, 121.5, 0.83, 1.19, 0.96, 5.66, 0.38, -0.10, 0.29],
        "V": [-2.69, -2.53, -1.29, 0, 5.9, 1.09, 0.13, 139.0, 1.06, 1.70, 0.50, 5.96, -0.75, -0.19, 0.03],
        "W": [-4.75, 3.65, 0.85, 1, 5.4, 0.88, 2.10, 228.2, 1.08, 1.37, 0.96, 5.89, -0.57, 0.31, 0.34],
        "Y": [-1.39, 2.32, 0.01, 1, 6.2, 0.33, 1.61, 197.0, 0.69, 1.47, 1.14, 5.66, -0.35, 0.40, -0.02],
    }

    def encode_single_signature(self, signature: str) -> List[float]:
        features = []
        for amino_acid_id in signature:
            properties = self.aa_params[amino_acid_id]
            features.extend(properties)
        return features

    def __call__(self, signatures: List[str]) -> np.ndarray:
        """_summary_

        Parameters
        ----------
        signatures : List[str]
            List of 34aa signatures.

        Returns
        -------
        np.ndarray of shape (n_samples, n_features)
            Paras features.
        """
        features = np.zeros((len(signatures), len(self.aa_params['-']) * 34))
        for i, signature in enumerate(signatures):
            if len(signature) != 34:
                raise ValueError("A-domain extended signature must be 34 amino acids long.")
            features[i] = self.encode_single_signature(signature)
        return features


class ParasWrapper(object):
    def __init__(self, 
                 model_dump: Path):
        """
        Parameters
        ----------
        model_dump : Path
            Path to joblib.dump'ed RandomForestClassifier model.
        """
        self.model = joblib.load(model_dump)
        self.featurizer = ParasFeaturizer()

    def predict_as_json(self, signatures: List[AA34], sort=False) -> List[Dict[str, Any]]:
        """_summary_

        Parameters
        ----------
        signatures : List[str]
            List of aa34 signatures
        sort : bool, optional
            Sort results by predicted probability, by default False

        Returns
        -------
        List[Dict[str, Any]]
            List of dicts following paras webserver output structure.
        """
        features = self.featurizer(signatures)
        predictions = self.model.predict_proba(features)
        class_labels = self.model.classes_

        def _form_prediction(substrate, prob):
            return {'probability': f'{prob:.3f}', 'substrate_name': substrate}

        results = []
        for i, signature in enumerate(signatures):
            cur_predictions = [_form_prediction(substrate, prob) 
                               for substrate, prob in zip(class_labels, predictions[i])]
            if sort:
                cur_predictions = sorted(cur_predictions, 
                                         key=lambda x: float(x['probability']), reverse=True)
            results.append({'domain_extended_signature': signature, 'predictions': cur_predictions})
        return results

    def predict(self, signature: AA34) -> Dict[PARAS_RESIDUE, Prob]:
        return {single_pred['substrate_name']: float(single_pred['probability'])
                for single_pred in self.predict_as_json([signature])[0]['predictions']}