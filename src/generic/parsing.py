from typing import Any, Callable, Optional

def read_int_pair(s: str) -> tuple[int, int]:
    # read a pair of ints from a string in the format "(a, b)"
    s = s.strip()
    if not (s.startswith('(') and s.endswith(')')):
        raise ValueError(f"Expected a pair of ints in the format '(a, b)', got '{s}'")
    s = s[1:-1]  # remove the parentheses
    parts = s.split(',')
    if len(parts) != 2:
        raise ValueError(f"Expected a pair of ints in the format '(a, b)', got '{s}'")
    try:
        a = int(parts[0].strip())
        b = int(parts[1].strip())
    except ValueError:
        raise ValueError(f"Expected a pair of ints in the format '(a, b)', got '{s}'")
    return a, b
