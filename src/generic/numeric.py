import math

def safe_log(x, base=math.e):
    """Returns log(x) or -inf if x <= 0."""
    if x <= 0:
        return float('-inf')
    return math.log(x, base)