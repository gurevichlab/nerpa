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

from src.antismash_parsing.genomic_context import ModuleGenomicContext, ModuleGenomicContextFeature
from src.antismash_parsing.bgc_variant_types import (
    BGC_ID,
    BGC_Variant_ID,
    BGC_Variant,
    GeneId,
)
from src.monomer_names_helper import NRP_Monomer
from src.general_type_aliases import (
    LogProb,
    Prob,
)
from src.monomer_names_helper import enum_representer, MonCode
import yaml


StateIdx = int
GenomicContext = Tuple[ModuleGenomicContextFeature, ...]

class HMM(NamedTuple):
    bgc_variant_id: BGC_Variant_ID
    transitions: List[List[Tuple[StateIdx, LogProb]]]  # u -> [(v, log_prob(u -> v))]
    emissions: List[List[LogProb]]  # u -> [log_prob(u -> emission)]
    module_start_states: List[StateIdx]
    module_match_states: List[StateIdx]

    def to_json(self):
        _emissions = [[emission if emission != LogProb(float('-inf')) else None
                      for emission in state_emissions]
                     for state_emissions in self.emissions]
        return {
            'bgc_variant_id': self.bgc_variant_id.to_dict(),
            'transitions': self.transitions,
            'emissions': _emissions,
            'module_start_states': self.module_start_states,
            'module_match_states': self.module_match_states
        }

# HMM_LOKC: emisson scores are Log-Odds,
# Unknown Chiralities are allowed,
# These HMMs are used for matching
# Transition and emission scores do NOT sum to 1
HMM_LOUC = NewType('HMM_LOUC', HMM)

# HMM_LPUC: emisson scores are Log-Probabilities,
# Unknown Chiralities are NOT allowed (the corresponding scores are set to -inf),
# These HMMs are used for p-value estimation
# Transition and emission scores sum to 1
HMM_LPKC = NewType('HMM_LPKC', HMM)

# HMM_LPUC: emisson scores are Log-Probabilities,
# Unknown Chiralities are allowed
# These HMMs are used for matching
# emission scores do NOT sum to 1
HMM_LPUC = NewType('HMM_LPUC', HMM)

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

    def __lt__(self, other):
        return self.name < other.name

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

