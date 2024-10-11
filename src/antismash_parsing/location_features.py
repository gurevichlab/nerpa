from typing import NewType, Tuple
from enum import Enum, auto

class ModuleLocFeature(Enum):
    START_OF_BGC = auto()
    END_OF_BGC = auto()
    START_OF_FRAGMENT = auto()
    END_OF_FRAGMENT = auto()
    START_OF_GENE = auto()
    END_OF_GENE = auto()
    PKS_UPSTREAM = auto()
    PKS_DOWNSTREAM = auto()


class GeneLocFeature(Enum):
    START_OF_BGC = auto()
    END_OF_BGC = auto()
    START_OF_FRAGMENT = auto()
    END_OF_FRAGMENT = auto()
    PKS_UPSTREAM = auto()
    PKS_DOWNSTREAM = auto()


class BGC_Fragment_Loc_Feature(Enum):
    START_OF_BGC = auto()
    END_OF_BGC = auto()
    PKS_UPSTREAM = auto()
    PKS_DOWNSTREAM = auto()

ModuleLocFeatures = Tuple[ModuleLocFeature, ...]
GeneLocFeatures = Tuple[GeneLocFeature, ...]
BGC_Fragment_Loc_Features = Tuple[BGC_Fragment_Loc_Feature, ...]