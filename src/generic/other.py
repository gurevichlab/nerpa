from typing import (
    Dict,
    List,
    NewType,
    Optional,
    Tuple,
    TypeVar
)
import re

T = TypeVar('T')

def get_score(scores: Dict[Tuple[T, ...], float], features: Tuple[T, ...]) -> float:
    if features in scores:
        return scores[features]
    features_set = set(features)
    return max(scores[x]
               for x in scores
               if set(x).issubset(features_set))


def parse_pretty_table(table: str) -> List[Dict[str, str]]:
    lines = table.splitlines()
    keys = lines[0].split()
    return [dict(zip(keys, line.split()))
            for line in lines[1:]]