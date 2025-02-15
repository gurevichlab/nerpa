from typing import List, Callable
from itertools import chain
import functools
import time


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


def make_optional(f: Callable) -> Callable:
    def optional_f(*args, **kwargs):
        # If any positional argument is `None`, return `None`
        if any(arg is None for arg in args):
            return None
        # Call the wrapped function
        return f(*args, **kwargs)
    return optional_f


def timing_decorator(func):
    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        start_time = time.time()  # Start timing
        result = func(self, *args, **kwargs)  # Execute the function
        elapsed_time = time.time() - start_time  # Calculate elapsed time

        # Use self.log.info if available; otherwise, fallback to print
        if hasattr(self, 'log') and hasattr(self.log, 'info'):
            self.log.info(f"Elapsed time: {elapsed_time:.4f} seconds")
        else:
            print(f"Elapsed time: {elapsed_time:.4f} seconds")

        return result

    return wrapper
