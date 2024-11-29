from typing import List, Tuple, NamedTuple, Dict, Optional
from dataclasses import dataclass
from enum import Enum, auto

from src.data_types import NRP_Monomer, BGC_Module, BGC_Variant


class HMM(NamedTuple):
    adj_list: List[List[Tuple[int, float]]]  # u -> [(v, log_prob(u -> v))]
    emission_log_probs: List[List[float]]  # u -> [log_prob(u -> emission)]


class DetailedHMMStateType(Enum):
    INITIAL = auto()

    SKIP_AT_BEGINNING_START = auto()
    SKIP_FRAGMENT_AT_BEGINNING = auto()

    MODULE_START = auto()
    MATCH = auto()
    INSERT_MONOMER = auto()

    SKIP_AT_END_START = auto()
    SKIP_FRAGMENT_END = auto()

    FINAL = auto()


@dataclass
class DetailedHMMState:
    state_type: DetailedHMMStateType
    emissions: Dict[NRP_Monomer, float]


class DetailedHMMEdgeType(Enum):
    START_SKIPPING_AT_START = auto()
    CONTINUE_SKIPPING_AT_START = auto()
    END_SKIPPING_AT_START = auto()

    MATCH = auto()
    NO_INSERTIONS = auto()
    START_INSERTING = auto()
    CONTINUE_INSERTING = auto()
    END_INSERTING = auto()

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

