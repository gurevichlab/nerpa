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
    StateIdx,
)
from src.hmm.hmm import HMM, HMM_LPUC
from src.matching.alignment_to_path_in_hmm import alignment_to_hmm_path
from src.matching.alignment_type import AlignmentLight, alignment_to_light_alignment
from src._not_used.p_values_estimation import PValueEstimator
from src.monomer_names_helper import Chirality
from src.rban_parsing.rban_monomer import rBAN_Monomer

from src.hmm.viterbi_algorithm import get_opt_path_with_score
from src.matching.hmm_to_alignment import hmm_path_to_alignment
from src.matching.alignment_type import Alignment
from src.hmm.hmm_scoring_helper import HMMHelper, ScoringScheme
from graphviz import Digraph
from pathlib import Path
from more_itertools import split_when


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
    _hmms: Dict[ScoringScheme, HMM] = None  # cache for HMMs with different emission weights and unk chirality settings
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
               scoring_scheme: ScoringScheme = ScoringScheme.NERPA) -> HMM:
        ST = DetailedHMMStateType  # for convenience

        if self._hmms is None:
            self._hmms = {}
        if scoring_scheme in self._hmms:
            return self._hmms[scoring_scheme]
        num_states = len(self.states)

        adj_list = [[(edge_to, edge_data.weight)
                     for edge_to, edge_data in self.transitions[edge_from].items()]
                    for edge_from in range(num_states)]

        def get_emissions(state_idx: StateIdx) -> Dict[NRP_Monomer, LogProb]:
            state = self.states[state_idx]
            match state.state_type:
                case (ST.MATCH
                      | ST.MATCH_POSSIBLE_ASSEMBLY_FINISH
                      | ST.MATCH_LAST_MODULE):
                    bgc_module_idx = self.state_idx_to_module_idx[state_idx]
                    bgc_module = self.bgc_variant.modules[bgc_module_idx]
                    return self.hmm_helper.get_match_emissions(
                        bgc_module=bgc_module,
                        scoring_scheme=scoring_scheme,
                        pks_domains_in_bgc=self.bgc_variant.has_pks_domains()
                    )
                case (ST.INSERT
                      | ST.INSERT_AT_START
                      | ST.INSERT_AT_END):
                    return self.hmm_helper.get_insert_emissions(
                        bgc_module=None,
                        scoring_scheme=scoring_scheme,
                        pks_domains_in_bgc=self.bgc_variant.has_pks_domains()
                    )
                case _:
                    return dict()
                                                               
        emission_scores = [
            [
                score
                for _mon, score in sorted(  # sort by monomer code to ensure consistent order of emissions across states
                        get_emissions(state_idx).items(),
                        key=lambda m_s: self.hmm_helper.monomer_names_helper.mon_to_int[m_s[0]]
                )
            ]
            for state_idx in range(len(self.states))
        ]

        module_names = [
            f'F{module.fragment_idx}:{module.gene_id}:{module.a_domain_idx}'
            for module in self.bgc_variant.modules
        ]

        self._hmms[scoring_scheme] = \
            HMM(
                bgc_variant_id=self.bgc_variant.bgc_variant_id,
                transitions=adj_list,
                emissions=emission_scores,
                module_start_states=self._module_idx_to_state_idx,
                module_match_states=self._module_idx_to_match_state_idx,
                state_types=[
                    state.state_type for state in self.states
                ],
                module_names=module_names
            )
        return self._hmms[scoring_scheme]

    def get_opt_path_with_emissions(self,
                                    start_state: StateIdx,
                                    finish_state: StateIdx,
                                    emitted_monomers: List[rBAN_Monomer]) -> List[Tuple[int, Optional[rBAN_Monomer]]]:
        monomer_codes = [self.hmm_helper.monomer_names_helper.mon_to_int[mon.to_base_mon()]
                         for mon in emitted_monomers]
        score, path = get_opt_path_with_score(hmm=self.to_hmm(),
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
             highlight_path: Optional[List[int]] = None,
             emission_names: Optional[List[str]] = None):
        self.to_hmm().draw(
            filename=filename,
            highlight_path=highlight_path,
            emission_names=emission_names
        )
