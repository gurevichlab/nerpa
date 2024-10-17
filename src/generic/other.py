from typing import (
    Dict,
    List,
    NewType,
    Optional,
    Tuple,
    TypeVar
)


T = TypeVar('T')

def get_score(scores: Dict[Tuple[T, ...], float], features: Tuple[T, ...]) -> float:
    if features in scores[features]:
        return scores[features]
    features_set = set(features)
    return max(scores[x]
               for x in scores
               if set(x).issubset(features_set))