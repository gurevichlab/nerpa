from __future__ import annotations
from typing import (
    List,
    Tuple,
    NamedTuple,
    Dict,
    Optional,
    TYPE_CHECKING,
    Union, NewType
)
from dataclasses import dataclass
from enum import Enum, auto

from src.antismash_parsing.genomic_context import ModuleGenomicContext
from src.data_types import (
    BGC_ID,
    BGC_Variant_ID,
    GeneId,
    NRP_Monomer,
    LogProb,
)
from src.monomer_names_helper import enum_representer, MonCode
import yaml


StateIdx = int

class HMM(NamedTuple):
    bgc_variant_id: BGC_Variant_ID
    transitions: List[List[Tuple[StateIdx, LogProb]]]  # u -> [(v, log_prob(u -> v))]
    emissions: List[List[LogProb]]  # u -> [log_prob(u -> emission)]
    module_start_states: List[StateIdx]
    module_match_states: List[StateIdx]

    def to_json(self):
        return {
            'bgc_variant_id': self.bgc_variant_id.to_dict(),
            'transitions': self.transitions,
            'emissions': self.emissions,
            'module_start_states': self.module_start_states,
            'module_match_states': self.module_match_states
        }

HMM_wo_unk_chir = NewType('HMM_wo_unk_chir', HMM)

class DetailedHMMStateType(Enum):
    INITIAL = auto()

    START_SKIP_MODULES_AT_START = auto()
    SKIP_MODULE_AT_START = auto()
    INSERT_AT_START = auto()

    MODULE_START = auto()
    MATCH = auto()
    INSERT = auto()
    CHOOSE_IF_ITERATE = auto()

    END_MATCHING = auto()
    INSERT_AT_END = auto()
    END_INSERTING_AT_END = auto()
    SKIP_MODULE_AT_END = auto()

    FINAL = auto()

    def is_emitting(self) -> bool:
        return self in {DetailedHMMStateType.MATCH,
                        DetailedHMMStateType.INSERT,
                        DetailedHMMStateType.INSERT_AT_START}


@dataclass
class DetailedHMMState:
    state_type: DetailedHMMStateType
    emissions: Dict[NRP_Monomer, LogProb]


class DetailedHMMEdgeType(Enum):
    START_INSERTING_AT_START = auto()
    INSERT_AT_START = auto()
    START_SKIP_MODULES_AT_START = auto()
    SKIP_MODULE_AT_START = auto()

    START_MATCHING = auto()

    MATCH = auto()
    NO_INSERTIONS = auto()
    START_INSERTING = auto()
    INSERT = auto()
    END_INSERTING = auto()

    ITERATE_MODULE = auto()
    ITERATE_GENE = auto()
    NO_ITERATION = auto()

    SKIP_MODULE = auto()
    SKIP_GENE = auto()

    END_MATCHING = auto()
    SKIP_MODULE_END_MATCHING = auto()
    START_INSERTING_AT_END = auto()
    INSERT_AT_END = auto()
    SKIP_MODULE_AT_END = auto()

yaml.add_representer(DetailedHMMEdgeType, enum_representer)

GenomicContext = ModuleGenomicContext

ModuleLevelEdgeKey = NewType('ModuleLevelEdgeKey', Tuple[GeneId, int])
GeneLevelEdgeKey = NewType('GeneLevelEdgeKey', GeneId)

# edge_key can be None for auxiliary edges that are not related to any BGC modules
EdgeKey = Union[ModuleLevelEdgeKey, GeneLevelEdgeKey, None]

class DetailedHMMEdge(NamedTuple):
    edge_type: DetailedHMMEdgeType
    weight: LogProb
    # edges have different weights, depending on the type and the context
    genomic_context: GenomicContext
    # edge_key is used in parameter estimation to not count the same edge multiple times
    # when it's used for different NRP compounds and the same BGC (but maybe different variants)
    edge_key: EdgeKey

