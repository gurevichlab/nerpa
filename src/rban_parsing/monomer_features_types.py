from typing import Tuple
from enum import Enum, auto


class MonomerFeature(Enum):
    START_OF_FRAGMENT = auto()
    END_OF_FRAGMENT = auto()
    UNKNOWN_RESIDUE = auto()
    PKS_HYBRID = auto()
    PKS = auto()


MonomerFeatures = Tuple[MonomerFeature, ...]
