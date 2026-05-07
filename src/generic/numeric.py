import math

def safe_log(x, base=math.e):
    """Returns log(x) or -inf if x <= 0."""
    if x <= 0:
        return float('-inf')
    return math.log(x, base)

def clamp(x, min_value, max_value):
    """Clamps x to the range [min_value, max_value]."""
    return max(min(x, max_value), min_value)
