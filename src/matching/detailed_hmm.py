from __future__ import annotations

from cmath import e
from copy import deepcopy
from typing import (
    ClassVar,
    Dict,
    List,
    NamedTuple,
    Tuple,
    Optional, Literal
)
from dataclasses import dataclass
from enum import Enum, auto

from src.data_types import (
    BGC_Module,
    BGC_Variant,
    BGC_Variant_ID,
    GeneId,
    LogProb,
    NRP_Monomer,
)
from src.matching.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    DetailedHMMStateType,
    DetailedHMMState,
    DetailedHMMEdge,
    HMM,
    StateIdx,
    GenomicContext,
    HMM_LOUC, HMM_LPKC, HMM_LPUC
)
from src.matching.alignment_to_path_in_hmm import alignment_to_hmm_path
from src.build_output.draw_hmm import draw_hmm
from src.matching.p_values_estimation import PValueEstimator
from src.monomer_names_helper import Chirality
from src.rban_parsing.rban_monomer import rBAN_Monomer

from src.matching.viterbi_algorithm import get_opt_path_with_score
from src.matching.bgc_to_hmm import bgc_variant_to_detailed_hmm
from src.matching.hmm_to_alignment import hmm_path_to_alignment
from src.matching.alignment_type import Alignment
from src.matching.hmm_checkpoints_heuristic import get_checkpoints
from src.matching.hmm_scoring_helper import HMMHelper
from src.matching.match_type import BGC_Variant_ID
from itertools import pairwise
from graphviz import Digraph
from pathlib import Path
from io import StringIO
from functools import cache


@dataclass
class DetailedHMM:
    states: List[DetailedHMMState]
    transitions: List[Dict[StateIdx, DetailedHMMEdge]]
    start_state_idx: StateIdx
    final_state_idx: StateIdx
    bgc_variant: BGC_Variant
    state_idx_to_module_idx: Dict[StateIdx, int]
    _module_idx_to_state_idx: List[StateIdx]  # points to MODULE_START state for each module. Used for building hmm from alignment
    _module_idx_to_match_state_idx: List[StateIdx]  # points to MATCH state for each module. Used for checkpoints heuristic

    hmm_helper: HMMHelper  # should be class variable but that would break parallelization for some reason
    _hmms: Dict[Tuple[Literal['LogProb', 'LogOdds'], bool], HMM] = None  # cache for HMMs with different emission weights and unk chirality settings
    _p_value_estimator: Optional[PValueEstimator] = None
    _score_vs_avg_nrp: Optional[LogProb] = None

    def score_vs_avg_bgc(self, nrp_monomers: List[rBAN_Monomer]) -> LogProb:
        return sum(self.hmm_helper.monomer_default_score(mon.to_base_mon())
                   for mon in nrp_monomers)

    def score_vs_avg_nrp(self) -> LogProb:
        '''
        Compute the null hypothesis score of the HMM, i.e.,
        the expected score of matching
        the BGC to an "average" NRP of the same length.
        '''
        if self._score_vs_avg_nrp is not None:
            return self._score_vs_avg_nrp

        total = 0
        for match_state in filter(lambda st: st.state_type == DetailedHMMStateType.MATCH,
                                  self.states):
            avg = 0
            for mon in self.hmm_helper.monomer_names_helper.proper_monomers():
                mon_default_freq = self.hmm_helper.monomer_names_helper.monomer_default_freq(mon)
                avg += mon_default_freq * match_state.emissions[mon]

            total += avg

        self._score_vs_avg_nrp = LogProb(total)
        return self._score_vs_avg_nrp

    @classmethod
    def from_bgc_variant(cls,
                         bgc_variant: BGC_Variant,
                         hmm_helper: HMMHelper) -> DetailedHMM:
        #print('from_bgc_variant', bgc_variant.genome_id)
        try:
            detailed_hmm = bgc_variant_to_detailed_hmm(DetailedHMM,
                                                       bgc_variant,
                                                       hmm_helper)
            return detailed_hmm
        except Exception as e:
            print(f'Failed to create DetailedHMM from BGC variant {bgc_variant.bgc_variant_id}')
            raise e

    # TODO: maybe rename unknown_chirality_allowed
    #  to improper_monomers_allowed or smth like that
    #  because it also affects PKS hybrids
    def to_hmm(self,
               emission_weights_type: Literal['LogProb', 'LogOdds']='LogProb',
               unknown_chirality_allowed: bool=False) -> HMM:
        if self._hmms is None:
            self._hmms = {}
        if (emission_weights_type, unknown_chirality_allowed) in self._hmms:
            return self._hmms[(emission_weights_type, unknown_chirality_allowed)]
        num_states = len(self.states)

        adj_list = [[(edge_to, edge_data.weight)
                     for edge_to, edge_data in self.transitions[edge_from].items()]
                    for edge_from in range(num_states)]

        def adjust_score(mon: NRP_Monomer, lp_score: LogProb) -> LogProb:
            if (not unknown_chirality_allowed and
                    (mon.chirality == Chirality.UNKNOWN or mon.is_pks_hybrid)):
                return LogProb(float('-inf'))
            match emission_weights_type:
                case 'LogProb': return lp_score
                case 'LogOdds': return lp_score - self.hmm_helper.monomer_default_score(mon)
                case _: raise ValueError(f'Invalid emission_weights_type: {emission_weights_type}. ')


        emission_scores = [[adjust_score(mon, state.emissions[mon])
                               for mon in sorted(state.emissions,
                                                 key=lambda m: self.hmm_helper.monomer_names_helper.mon_to_int[m])]  # TODO: define int(m) or use monomer_names_helper
                              for state in self.states]
        self._hmms[(emission_weights_type, unknown_chirality_allowed)] = \
            HMM(transitions=adj_list,
                emissions=emission_scores,
                module_start_states=self._module_idx_to_state_idx,
                module_match_states=self._module_idx_to_match_state_idx,
                bgc_variant_id=self.bgc_variant.bgc_variant_id)
        return self._hmms[(emission_weights_type, unknown_chirality_allowed)]

    def to_hmm_lpuc(self) -> HMM_LPUC:
        return HMM_LPUC(self.to_hmm(emission_weights_type='LogProb',
                                    unknown_chirality_allowed=True))

    def get_p_value(self, lo_score: LogProb, lp_score: LogProb) -> float:
        return self._p_value_estimator(lo_score, lp_score)

    def get_opt_path_with_emissions(self,
                                    start_state: StateIdx,
                                    finish_state: StateIdx,
                                    emitted_monomers: List[rBAN_Monomer]) -> List[Tuple[int, Optional[rBAN_Monomer]]]:
        monomer_codes = [self.hmm_helper.monomer_names_helper.mon_to_int[mon.to_base_mon()]
                         for mon in emitted_monomers]
        score, path = get_opt_path_with_score(hmm=self.to_hmm(emission_weights_type='LogProb',
                                                              unknown_chirality_allowed=True),
                                              observed_sequence=monomer_codes,
                                              start_state=start_state,
                                              finish_state=finish_state)

        emitted_monomers_iter = iter(emitted_monomers)
        path_with_emissions = []
        for state_idx in path:
            if self.states[state_idx].state_type.is_emitting():
                path_with_emissions.append((state_idx, next(emitted_monomers_iter)))
            else:
                path_with_emissions.append((state_idx, None))
        return path_with_emissions

    def path_to_alignment(self,
                          path: List[int],
                          nrp_monomers: List[rBAN_Monomer]) -> Alignment:
        #self.draw(Path(f'{self.bgc_variant.bgc_variant_id.bgc_id.genome_id}.png'), path)
        return hmm_path_to_alignment(self, path, nrp_monomers)

    def alignment_to_path_with_emisions(self, alignment: Alignment) -> List[Tuple[int, Optional[NRP_Monomer]]]:
        return alignment_to_hmm_path(self, alignment)

    def draw(self,
             filename: Path,
             highlight_path: Optional[List[int]] = None) -> Digraph:
        return draw_hmm(self, filename, highlight_path)


# captures all context relevant for parameter training
class EdgeFingerprint(NamedTuple):
    edge_type: DetailedHMMEdgeType
    edge_context: Optional[GenomicContext] = None

'''
def edge_fingerprint(hmm: DetailedHMM,
                     path: List[Tuple[int, NRP_Monomer]],
                     edge_idx: int) -> EdgeFingerprint:
    u, v = path[edge_idx][0], path[edge_idx + 1][0]
    edge = hmm.transitions[u][v]

    # find corresponding bgc module
    match edge.edge_type:
        # genomic info not relevant
        case (
        # before matching
            DetailedHMMEdgeType.START_INSERTING_AT_START |
            DetailedHMMEdgeType.INSERT_AT_START |
            DetailedHMMEdgeType.START_SKIP_MODULES_AT_START |
            DetailedHMMEdgeType.START_SKIP_GENES_AT_START |
            DetailedHMMEdgeType.START_SKIP_FRAGMENTS_AT_START |
        # at matching
            DetailedHMMEdgeType.INSERT |
            DetailedHMMEdgeType.END_INSERTING
        ):
            return EdgeFingerprint(edge_type=edge.edge_type)
        # only fragment context relevant
        case (
            DetailedHMMEdgeType.SKIP_FRAGMENT_AT_START |
            DetailedHMMEdgeType.SKIP_FRAGMENT |
            DetailedHMMEdgeType.SKIP_FRAGMENT_AT_END
        ):
            return EdgeFingerprint(edge_type=edge.edge_type,
                                   bgc_fragment_context=hmm.bgc_variant.fragments[hmm.state_idx_to_module_idx[u]].loc)
'''

