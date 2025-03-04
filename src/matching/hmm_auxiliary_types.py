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
from src.data_types import NRP_Monomer, LogProb
from src.monomer_names_helper import enum_representer, MonCode
from src.matching.match_type import Match_BGC_Variant_Info
import yaml


StateIdx = int

class HMM(NamedTuple):
    bgc_info: Match_BGC_Variant_Info
    transitions: List[List[Tuple[StateIdx, LogProb]]]  # u -> [(v, log_prob(u -> v))]
    emissions: List[List[LogProb]]  # u -> [log_prob(u -> emission)]
    module_start_states: List[StateIdx]
    module_match_states: List[StateIdx]

    # TODO: fix the discrepancy in naming: adj_list vs transitions, emission_log_probs vs emissions
    def to_json(self):
        return {
            'bgc_info': self.bgc_info._asdict(),
            'transitions': self.transitions,
            'emissions': self.emissions,
            'module_start_states': self.module_start_states,
            'module_match_states': self.module_match_states
        }


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

    def is_emitting(self) -> bool:
        return self in {DetailedHMMStateType.MATCH,
                        DetailedHMMStateType.INSERT,
                        DetailedHMMStateType.INSERT_AT_START}


@dataclass
class DetailedHMMState:
    state_type: DetailedHMMStateType
    emissions: Dict[NRP_Monomer, LogProb]


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
    weight: LogProb
    # edges have different weights, depending on the type and the context
    genomic_context: Optional[GenomicContext]
    # edge_key is used in parameter estimation to not count the same edge multiple times when it's used for different NRP compounds and the same BCG
    edge_key: Optional[tuple]

