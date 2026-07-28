from __future__ import annotations
from typing import (
    Tuple,
    NamedTuple,
    Dict,
    Optional,
)
from dataclasses import dataclass
from enum import Enum, auto

from src.antismash_parsing.genomic_context import ModuleGenomicContextFeature
from src.monomer_names_helper import NRP_Monomer
from src.general_type_aliases import (
    LogProb,
)
from src.monomer_names_helper import enum_representer
import yaml


StateIdx = int
GenomicContext = Tuple[ModuleGenomicContextFeature, ...]


class DetailedHMMStateType(Enum):
    INITIAL = auto()

    SKIPPING_MODULES_AT_START_STARTED = auto()
    SKIPPING_MODULES_AT_START = auto()
    INSERT_AT_START = auto()
    SKIPPING_MODULES_AT_START_FINISHED = auto()
    CHOOSE_IF_START_MATCHING = auto()
    SKIP_MODULE_AT_START = auto()  # a phantom state -- I pretend there's a chain of these states when computing skip costs

    MODULE_SUBGRAPH_ROOT = auto()
    MATCH = auto()
    MATCH_POSSIBLE_ASSEMBLY_FINISH = auto()
    MATCH_LAST_MODULE = auto()
    CHOOSE_IF_FINISH_MATCHING = auto()
    MATCHING_CONTINUED = auto()
    INSERT = auto()
    CHOOSE_IF_ITERATE_MODULE = auto()
    CHOOSE_IF_ITERATE_GENE = auto()

    MATCHING_FINISHED = auto()
    INSERT_AT_END = auto()
    INSERTING_AT_END_FINISHED = auto()
    SKIPPING_MODULES_AT_END = auto()
    SKIP_MODULE_AT_END = auto()  # a phantom state -- I pretend there's a chain of these states when computing skip costs

    CHOOSE_IF_ITERATE_BGC = auto()

    FINAL = auto()


@dataclass
class DetailedHMMState:
    state_type: DetailedHMMStateType
    emissions: Dict[NRP_Monomer, LogProb] = None
    related_module_idx: Optional[int] = None

    def __init__(self,
                 state_type: DetailedHMMStateType,
                 emissions: Optional[Dict[NRP_Monomer, LogProb]] = None,
                 related_module_idx: Optional[int] = None):
        self.state_type = state_type
        self.emissions = emissions if emissions is not None else {}
        self.related_module_idx = related_module_idx

    def is_emitting(self) -> bool:
        return self.emissions is not None and len(self.emissions) > 0


class DetailedHMMEdgeType(Enum):
    INITIAL = auto()
    START_INSERTING_AT_START = auto()
    INSERT_AT_START = auto()
    START_SKIPPING_MODULES_AT_START = auto()
    SKIP_MODULE_AT_START = auto()
    SKIP_UNTIL_NEXT_TENTATIVE_ASSEMBLY_LINE_START = auto()

    START_MATCHING = auto()

    MATCH = auto()
    NO_INSERTIONS = auto()
    START_INSERTING = auto()
    INSERT = auto()
    FINISH_INSERTING = auto()

    ITERATE_MODULE = auto()
    ITERATE_GENE = auto()
    NO_ITERATION = auto()
    CHOOSE_IF_FINISH_MATCHING = auto()
    FINISH_MATCHING = auto()
    CONTINUE_MATCHING = auto()

    SKIP_MODULE = auto()
    SKIP_GENE = auto()

    END_MATCHING = auto()
    START_SKIPPING_MODULES_AT_END = auto()
    SKIP_MODULE_END_MATCHING = auto()
    START_INSERTING_AT_END = auto()
    INSERT_AT_END = auto()
    SKIP_MODULES_AT_END = auto()
    SKIP_MODULE_AT_END = auto()

    ITERATE_BGC = auto()

    def __lt__(self, other):
        return self.name < other.name

    def is_insert_like(self) -> bool:
        return self in {
            DetailedHMMEdgeType.START_INSERTING_AT_START,
            DetailedHMMEdgeType.INSERT_AT_START,
            DetailedHMMEdgeType.START_INSERTING,
            DetailedHMMEdgeType.INSERT,
            DetailedHMMEdgeType.START_INSERTING_AT_END,
            DetailedHMMEdgeType.INSERT_AT_END,
        }

yaml.add_representer(DetailedHMMEdgeType, enum_representer)

#GenomicContext = ModuleGenomicContext

#ModuleLevelEdgeKey = NewType('ModuleLevelEdgeKey', Tuple[GeneId, int])
#GeneLevelEdgeKey = NewType('GeneLevelEdgeKey', GeneId)

# edge_key can be None for auxiliary edges that are not related to any BGC modules
# TODO: I think it's not used anywhere, can be removed
#EdgeKey = Union[ModuleLevelEdgeKey, GeneLevelEdgeKey, None]

class DetailedHMMEdge(NamedTuple):
    edge_type: DetailedHMMEdgeType
    weight: LogProb = 1  # default positive value for invalid log probability
    # edges have different weights, depending on the type and the context
    #genomic_context: GenomicContext
    # edge_key is used in parameter estimation to not count the same edge multiple times
    # when it's used for different NRP compounds and the same BGC (but maybe different variants)
    #edge_key: EdgeKey

