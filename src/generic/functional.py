from typing import List, Callable, Optional
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


def timing_decorator(task_name: Optional[str] = None):
    """Decorator that measures execution time and logs it."""

    # Handle both @timing_decorator and @timing_decorator("Some Task")
    if callable(task_name):  # Case: @timing_decorator used without parentheses
        func = task_name
        task_name = None  # No explicit task name was provided
        return timing_decorator()(func)  # Reapply with proper structure

    def _timing_decorator(func: Callable):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            start_time = time.time()
            result = func(*args, **kwargs)
            elapsed_time = time.time() - start_time

            # Create log message
            message = f"{task_name} took {elapsed_time:.4f} seconds" if task_name else f"Elapsed time: {elapsed_time:.4f} seconds"

            # Check if the first argument is an instance with `log.info`
            if args and hasattr(args[0], 'log') and hasattr(args[0].log, 'info'):
                args[0].log.info(message)
            else:
                print(message)

            return result

        return wrapper

    return _timing_decorator
