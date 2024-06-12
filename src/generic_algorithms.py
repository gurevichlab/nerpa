from typing import List, Callable
from itertools import chain


def hamming_distance(s1: str, s2: str) -> int:
    if len(s1) != len(s2):
        raise ValueError("Strings must be of equal length")
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def list_monad_compose(*functions: Callable) -> Callable[[List], List]:
    def composed(lst: List) -> List:
        for f in functions:
            lst = list(chain(*(map(f, lst))))
        return lst
    return composed
