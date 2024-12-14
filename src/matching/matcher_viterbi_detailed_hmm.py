from __future__ import annotations
from typing import List, Tuple, NamedTuple, Dict, Optional
from dataclasses import dataclass
from enum import Enum, auto

from src.data_types import NRP_Monomer, BGC_Module, BGC_Variant, GeneId, LogProb
from src.matching.matcher_viterbi_types import (
    DetailedHMMEdgeType,
    DetailedHMMStateType,
    DetailedHMMState,
    DetailedHMMEdge,
    HMM
)
from src.matching.matcher_viterbi_algorithm import get_opt_path_with_emissions
from src.matching.alignment_to_path_in_hmm import alignment_to_hmm_path
from src.build_output.draw_hmm import draw_hmm
from collections import defaultdict

from src.matching.bgc_to_hmm import bgc_variant_to_detailed_hmm
from src.matching.hmm_to_alignment import hmm_path_to_alignment
from src.matching.matching_types_alignment import Alignment
from src.matching.scoring_helper import ScoringHelper
from itertools import pairwise
from graphviz import Digraph
from pathlib import Path


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
                               for mon in sorted(state.emissions, key=lambda m: int(m))]  # TODO: define int(m) or use monomer_names_helper
                              for state in self.states]
        return HMM(adj_list=adj_list,
                   emission_log_probs=emission_log_probs)

    def get_opt_path_with_emissions(self,
                                    start: int,
                                    finish: int,
                                    emitted_monomers: List[NRP_Monomer]) -> List[Tuple[int, Optional[NRP_Monomer]]]:
        score, path_with_emissions = get_opt_path_with_emissions(self.to_hmm(),
                                                                 start, finish,
                                                                 [int(mon) for mon in emitted_monomers])
        return [(u, NRP_Monomer.from_int(mon_int))
                for u, mon_int in path_with_emissions]

    def path_to_alignment(self,
                          path: List[int],
                          nrp_monomers: List[NRP_Monomer],
                          scoring_helper: ScoringHelper) -> Alignment:
        return hmm_path_to_alignment(self, path, nrp_monomers, scoring_helper)

    def alignment_to_path_with_emisions(self, alignment: Alignment) -> List[Tuple[int, Optional[NRP_Monomer]]]:
        return alignment_to_hmm_path(self, alignment)

    def draw(self, filename: Path) -> Digraph:
        return draw_hmm(self, filename)

