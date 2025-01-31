from __future__ import annotations
from typing import (
    List,
    Tuple,
    NamedTuple,
    Dict,
    Optional,
    Union
)
from dataclasses import dataclass
from enum import Enum, auto

from src.antismash_parsing.location_features import (
    ModuleLocFeatures,
    GeneLocFeatures,
    BGC_Fragment_Loc_Features
)
from src.data_types import NRP_Monomer
from src.monomer_names_helper import enum_representer
import yaml

class HMM(NamedTuple):
    adj_list: List[List[Tuple[int, float]]]  # u -> [(v, log_prob(u -> v))]
    emission_log_probs: List[List[float]]  # u -> [log_prob(u -> emission)]


class DetailedHMMStateType(Enum):
    INITIAL = auto()

    SKIP_MODULE_AT_START = auto()
    SKIP_GENE_AT_START = auto()
    SKIP_FRAGMENT_AT_START = auto()
    INSERT_AT_START = auto()

    MODULE_START = auto()
    MATCH = auto()
    INSERT = auto()

    SKIP_AT_END_START = auto()
    SKIP_FRAGMENT_END = auto()

    FINAL = auto()


@dataclass
class DetailedHMMState:
    state_type: DetailedHMMStateType
    emissions: Dict[NRP_Monomer, float]


# TODO: I don't like these "CONTINUE" states,
#  maybe refactor to have just INSERT, etc. To do this I need to revisit the HMM generation
#  -- so that first goes the emission and then the transition
class DetailedHMMEdgeType(Enum):
    START_INSERTING_AT_START = auto()
    INSERT_AT_START = auto()
    START_SKIP_MODULES_AT_START = auto()
    START_SKIP_GENES_AT_START = auto()
    START_SKIP_FRAGMENTS_AT_START = auto()
    SKIP_FRAGMENT_AT_START = auto()  # for skipping genes and modules at start, regular SKIP_GENE/SKIP_MODULE are used

    START_MATCHING = auto()

    MATCH = auto()
    NO_INSERTIONS = auto()
    START_INSERTING = auto()
    INSERT = auto()
    END_INSERTING = auto()

    ITERATE_MODULE = auto()
    ITERATE_GENE = auto()

    SKIP_MODULE = auto()
    SKIP_GENE = auto()
    SKIP_FRAGMENT = auto()

    START_SKIPPING_AT_END = auto()
    SKIP_FRAGMENT_AT_END = auto()

yaml.add_representer(DetailedHMMEdgeType, enum_representer)

GenomicContext = Union[ModuleLocFeatures, GeneLocFeatures, BGC_Fragment_Loc_Features]
EdgeKey = tuple

class DetailedHMMEdge(NamedTuple):
    edge_type: DetailedHMMEdgeType
    log_prob: float
    # edges have different weights, depending on the type and the context
    genomic_context: Optional[GenomicContext]
    # edge_key is used in parameter estimation to not count the same edge multiple times when it's used for different NRP compounds and the same BCG
    edge_key: Optional[tuple]

