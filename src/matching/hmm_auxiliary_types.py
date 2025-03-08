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

from src.antismash_parsing.genomic_context import (
    ModuleGenomicContext,
    GeneGenomicContext,
    FragmentGenomicContext
)
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
    bgc_info: BGC_Variant_ID
    transitions: List[List[Tuple[StateIdx, LogProb]]]  # u -> [(v, log_prob(u -> v))]
    emissions: List[List[LogProb]]  # u -> [log_prob(u -> emission)]
    module_start_states: List[StateIdx]
    module_match_states: List[StateIdx]

    def to_json(self):
        return {
            'bgc_info': self.bgc_info.to_dict(),
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

GenomicContext = Union[ModuleGenomicContext, GeneGenomicContext, FragmentGenomicContext]


ModuleLevelEdgeKey = NewType('ModuleLevelEdgeKey', Tuple[GeneId, int])
GeneLevelEdgeKey = NewType('GeneLevelEdgeKey', GeneId)
FragmentLevelEdgeKey = NewType('FragmentLevelEdgeKey', Tuple[GeneId, ...])

# edge_key can be None for auxiliary edges that are not related to any BGC modules
EdgeKey = Union[ModuleLevelEdgeKey, GeneLevelEdgeKey, FragmentLevelEdgeKey, None]

class DetailedHMMEdge(NamedTuple):
    edge_type: DetailedHMMEdgeType
    weight: LogProb
    # edges have different weights, depending on the type and the context
    genomic_context: Optional[GenomicContext]
    # edge_key is used in parameter estimation to not count the same edge multiple times
    # when it's used for different NRP compounds and the same BGC (but maybe different variants)
    edge_key: EdgeKey

