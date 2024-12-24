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
from src.data_types import NRP_Monomer, BGC_Module, BGC_Variant
from src.matching.bgc_to_hmm import bgc_variant_to_detailed_hmm
from itertools import pairwise
from src.matching.matching_types_alignment_step import (
    AlignmentStep,
    AlignmentStep_BGC_Module_Info,
    AlignmentStep_NRP_Monomer_Info
)
from src.matching.matching_types_alignment import Alignment
from src.matching.scoring_helper import ScoringHelper


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


GenomicContext = Union[ModuleLocFeatures, GeneLocFeatures, BGC_Fragment_Loc_Features]
EdgeKey = tuple

class DetailedHMMEdge(NamedTuple):
    edge_type: DetailedHMMEdgeType
    log_prob: float
    # edges have different weights, depending on the type and the context
    genomic_context: Optional[GenomicContext]
    # edge_key is used in parameter estimation to not count the same edge multiple times when it's used for different NRP compounds and the same BCG
    edge_key: Optional[tuple]


'''
@dataclass
class DetailedHMM:
    states: List[DetailedHMMState]
    adj_list: List[List[DetailedHMMEdge]]
    start_state_idx: int
    final_state_idx: int
    bgc_variant: BGC_Variant
    state_idx_to_module_idx: Dict[int, int]
    _module_idx_to_state_idx: List[int]  # points to MODULE_START state for each module. Used for building hmm from alignment

    @classmethod
    def from_bgc_variant(cls, bgc_variant: BGC_Variant) -> DetailedHMM:
        return bgc_variant_to_detailed_hmm(DetailedHMM, bgc_variant)

    def to_hmm(self) -> HMM:
        num_states = len(self.states)
        adj_list = [[(edge.to, edge.log_prob) for edge in self.adj_list[u]]
                    for u in range(num_states)]

        emission_log_probs = [[state.emissions[mon]
                               for mon in sorted(state.emissions, key=lambda m: int(m))]
                              for state in self.states]
        return HMM(adj_list=adj_list,
                   emission_log_probs=emission_log_probs)

    def path_to_alignment(self,
                          path: List[int],
                          nrp_monomers: List[NRP_Monomer],
                          scoring_helper: ScoringHelper) -> Alignment:
        # The idea is that first the transition is made and then the emission -- in the new state
        # TODO: maybe refactor to do first emission and then transition to make it more similar to the standard Viterbi algorithm
        alignment = []
        mon_idx = 0
        for edge_from, edge_to in pairwise(path):
            edge = next(edge for edge in self.adj_list[edge_from] if edge.to == edge_to)  # maybe use dicts for speedup?

            # these edge types are handled separately because they involve multiple modules
            # and hence yield multiple alignment steps
            if edge.edge_type in (DetailedHMMEdgeType.SKIP_GENE,
                                  DetailedHMMEdgeType.SKIP_FRAGMENT,
                                  DetailedHMMEdgeType.SKIP_FRAGMENT_AT_START,
                                  DetailedHMMEdgeType.SKIP_FRAGMENT_AT_END):
                fst_module_idx = self.state_idx_to_module_idx[edge_from]
                last_module_idx = self.state_idx_to_module_idx[edge_to]
                alignment.extend(
                    AlignmentStep(bgc_module_info=AlignmentStep_BGC_Module_Info.from_bgc_module(self.bgc_variant.modules[module_idx]),
                                  nrp_monomer_info=None,
                                  score=edge.log_prob,
                                  match_detailed_score=None,
                                  step_type=edge.edge_type)
                    for module_idx in range(fst_module_idx, last_module_idx + 1))
                continue

            bgc_module = None
            nrp_monomer = None
            match_detailed_score = None
            transition_score = edge.log_prob
            emission_score = 0

            match edge.edge_type:
                # with emissions
                case DetailedHMMEdgeType.MATCH:
                    bgc_module = self.bgc_variant.modules[self.state_idx_to_module_idx[edge_from]]
                    nrp_monomer = nrp_monomers[mon_idx]
                    match_detailed_score = scoring_helper.match_detailed_score(bgc_module, nrp_monomers[mon_idx])
                    emission_score = sum(match_detailed_score)
                    mon_idx += 1
                case DetailedHMMEdgeType.START_INSERTING | DetailedHMMEdgeType.INSERT:
                    bgc_module = self.bgc_variant.modules[self.state_idx_to_module_idx[edge_from]]
                    nrp_monomer = nrp_monomers[mon_idx]
                    # in the current setup, emission_score is always 0 for insertions, as if the emission is always the same
                    mon_idx += 1
                case DetailedHMMEdgeType.START_INSERTING_AT_START | DetailedHMMEdgeType.INSERT_AT_START:
                    next_module_state_idx = next(next_edge.to
                                                 for next_edge in self.adj_list[edge_to]
                                                 if self.states[
                                                     next_edge.to].state_type == DetailedHMMStateType.MODULE_START)
                    bgc_module = self.bgc_variant.modules[next_module_state_idx]
                    # in the current setup, emission_score is always 0 for insertions, as if the emission is always the same
                    nrp_monomer = nrp_monomers[mon_idx]
                    mon_idx += 1
                # skips
                case DetailedHMMEdgeType.SKIP_MODULE | DetailedHMMEdgeType.SKIP_GENE | DetailedHMMEdgeType.SKIP_FRAGMENT:
                    bgc_module = self.bgc_variant.modules[self.state_idx_to_module_idx[edge_from]]

            score = transition_score + emission_score
            bgc_module_info = AlignmentStep_BGC_Module_Info.from_bgc_module(bgc_module) \
                if bgc_module is not None else None
            nrp_monomer_info = AlignmentStep_NRP_Monomer_Info.from_nrp_monomer(nrp_monomer) \
                if nrp_monomer is not None else None
            alignment.append(AlignmentStep(bgc_module_info=bgc_module_info,
                                           nrp_monomer_info=nrp_monomer_info,
                                           score=score,
                                           match_detailed_score=match_detailed_score,
                                           step_type=edge.edge_type))

        return alignment


'''