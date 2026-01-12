from __future__ import annotations

from typing import (
    Dict,
    List,
    Tuple,
    Optional, Literal
)
from dataclasses import dataclass

from src.antismash_parsing.bgc_variant_types import (
    BGC_Variant,
)
from src.hmm.hmm_constructor.hmm_constructor import HMM_Constructor
from src.monomer_names_helper import NRP_Monomer, PKS_RESIDUE, NOT_NRPS_RESIDUE
from src.general_type_aliases import (
    LogProb,
)
from src.hmm.hmm_auxiliary_types import (
    DetailedHMMStateType,
    DetailedHMMState,
    DetailedHMMEdge,
    HMM,
    StateIdx,
    HMM_LPUC
)
from src.matching.alignment_to_path_in_hmm import alignment_to_hmm_path
from src.matching.alignment_type import AlignmentLight, alignment_to_light_alignment
from src.build_output.draw_hmm import draw_hmm
from src._not_used.p_values_estimation import PValueEstimator
from src.monomer_names_helper import Chirality
from src.rban_parsing.rban_monomer import rBAN_Monomer

from src.hmm.viterbi_algorithm import get_opt_path_with_score
from src.matching.hmm_to_alignment import hmm_path_to_alignment
from src.matching.alignment_type import Alignment
from src.hmm.hmm_scoring_helper import HMMHelper
from graphviz import Digraph
from pathlib import Path


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

    def score_vs_avg_bgc(self,
                         nrp_monomers: List[rBAN_Monomer],
                         specificities_source: Literal['PARAS average', 'Norine frequencies'] = 'PARAS average') -> LogProb:
        def paras_avg_mon_score(_mon: rBAN_Monomer) -> LogProb:
            mon = _mon.to_base_mon()
            if mon.residue in (PKS_RESIDUE, NOT_NRPS_RESIDUE):
                insert_scores = self.hmm_helper.get_insert_emissions(bgc_module=None)
                return insert_scores[mon]
            else:
                return self.hmm_helper.monomer_default_score(mon)

        def norine_mon_freq_score(mon: rBAN_Monomer) -> LogProb:
            return self.hmm_helper.monomer_names_helper.monomer_default_log_freq(mon.to_base_mon())

        if specificities_source == 'PARAS average':
            mon_score = paras_avg_mon_score
        elif specificities_source == 'Norine frequencies':
            mon_score = norine_mon_freq_score
        else:
            raise ValueError(f'Invalid specificities_source: {specificities_source}')

        return sum(mon_score(mon) for mon in nrp_monomers)

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
        #print('Constructing HMM from bgc variant', bgc_variant.bgc_variant_id.bgc_id.to_str_short())
        # detailed_hmm = bgc_variant_to_detailed_hmm(DetailedHMM,
        #                                            bgc_variant,
        #                                            hmm_helper)
        detailed_hmm = HMM_Constructor(bgc_variant=bgc_variant,
                                       hmm_helper=hmm_helper).build_hmm(DetailedHMM)
        #detailed_hmm.draw(Path(f'detailed_hmm.png'))

        return detailed_hmm

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
            if self.states[state_idx].emissions:
                path_with_emissions.append((state_idx, next(emitted_monomers_iter)))
            else:
                path_with_emissions.append((state_idx, None))
        return path_with_emissions

    def path_to_alignment(self,
                          path: List[int],
                          nrp_monomers: List[rBAN_Monomer]) -> Alignment:
        #self.draw(Path(f'{self.bgc_variant.bgc_variant_id.bgc_id.genome_id}.png'), path)
        return hmm_path_to_alignment(self, path, nrp_monomers)

    # have to specify alignment type because Python can't infer types properly
    def alignment_to_path_with_emisions(self,
                                        alignment: Alignment | AlignmentLight,
                                        alignment_type: Literal['full', 'light'] = 'full') \
            -> List[Tuple[int, Optional[NRP_Monomer]]]:
        _alignment = (
            alignment_to_light_alignment(alignment)
            if alignment_type == 'full'
            else alignment
        )
        return alignment_to_hmm_path(self, _alignment)

    def draw(self,
             filename: Path,
             highlight_path: Optional[List[int]] = None) -> Digraph:
        return draw_hmm(self, filename, highlight_path)


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

