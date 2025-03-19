from typing import (
    Dict,
    List,
    NewType,
    Optional,
    Tuple,
    TypeVar
)
import re
from src.antismash_parsing.location_features import ModuleLocFeature

T = TypeVar('T')

def get_score(scores: Dict[Tuple[T, ...], float], features: Tuple[T, ...]) -> float:
    if features in scores:
        return scores[features]
    features_set = set(features)
    if ModuleLocFeature.END_OF_BGC in features_set:
        pass
    result = max(scores[x]
                 for x in scores
                 if set(x).issubset(features_set))
    return result


def parse_pretty_table(table: str) -> List[Dict[str, str]]:
    lines = table.splitlines()
    keys = lines[0].split()
    return [dict(zip(keys, line.split()))
            for line in lines[1:]]