from typing import (
    Dict,
    List,
    NamedTuple,
    NewType,
    Optional,
    Tuple,
    TypeVar
)
import re
from src.antismash_parsing.genomic_context import ModuleGenomicContextFeature

T = TypeVar('T')

def get_score(scores: Dict[Tuple[T, ...], float], features: Tuple[T, ...]) -> float:
    if features in scores:
        return scores[features]
    features_set = set(features)
    if ModuleGenomicContextFeature.END_OF_BGC in features_set:
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


def json_round_floats(obj, ndigits=3):
    '''
    Recursively round all floats in a JSON-like structure (dicts, lists, floats).
    '''
    if isinstance(obj, float):
        return round(obj, ndigits)
    elif isinstance(obj, dict):
        return {k: json_round_floats(v, ndigits) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [json_round_floats(v, ndigits) for v in obj]
    elif isinstance(obj, tuple):
        return tuple(json_round_floats(v, ndigits) for v in obj)
    else:
        return obj
