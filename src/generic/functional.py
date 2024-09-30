from typing import List, Callable
from itertools import chain


def list_monad_compose(*functions: Callable) -> Callable[[List], List]:
    def composed(lst: List) -> List:
        for f in functions:
            lst = list(chain(*(map(f, lst))))
        return lst
    return composed


class CachedByKey:
    f: Callable
    key: Callable
    cache: dict

    def __init__(self, f: Callable, key: Callable):
        self.f = f
        self.key = key
        self.cache = {}

    def __call__(self, *args, **kwargs):
        k = self.key(*args, **kwargs)
        if k not in self.cache:
            self.cache[k] = self.f(*args, **kwargs)
        return self.cache[k]


def cached_by_key(key: Callable):
    def decorator(f: Callable):
        return CachedByKey(f, key)
    return decorator


def compose(*funs):
    if len(funs) == 1:
        return funs[0]
    else:
        return funs[0](compose(*funs[1:]))
