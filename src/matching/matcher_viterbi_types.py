from typing import List, Tuple, NamedTuple, Dict, Optional
from dataclasses import dataclass
from enum import Enum, auto

from src.data_types import NRP_Monomer, BGC_Module, BGC_Variant


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


class DetailedHMMEdgeType(Enum):
    START_INSERTING_AT_START = auto()
    START_SKIP_MODULES_AT_START = auto()
    START_SKIP_GENES_AT_START = auto()
    START_SKIP_FRAGMENTS_AT_START = auto()
    CONTINUE_SKIP_FRAGMENTS_AT_START = auto()

    START_MATCHING = auto()

    MATCH = auto()
    NO_INSERTIONS = auto()
    START_INSERTING = auto()
    CONTINUE_INSERTING = auto()
    END_INSERTING = auto()

    ITERATE_MODULE = auto()
    ITERATE_GENE = auto()

    SKIP_MODULE = auto()
    SKIP_GENE = auto()
    SKIP_FRAGMENT = auto()

    START_SKIPPING_AT_END = auto()
    CONTINUE_SKIPPING_AT_END = auto()


class DetailedHMMEdge(NamedTuple):
    edge_type: DetailedHMMEdgeType
    log_prob: float
    to: int


@dataclass
class DetailedHMM:
    states: List[DetailedHMMState]
    adj_list: List[List[DetailedHMMEdge]]
    bgc_variant: BGC_Variant
    state_idx_to_module_idx: Dict[int, int]

