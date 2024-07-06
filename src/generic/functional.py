from typing import List, Callable
from itertools import chain


def list_monad_compose(*functions: Callable) -> Callable[[List], List]:
    def composed(lst: List) -> List:
        for f in functions:
            lst = list(chain(*(map(f, lst))))
        return lst
    return composed