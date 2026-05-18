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
from src.hmm.hmm_auxiliary_types import (
    DetailedHMMStateType,
    StateIdx,
)
from src.monomer_names_helper import NRP_Monomer
from src.general_type_aliases import (
    LogProb,
    Prob,
)
from src.monomer_names_helper import enum_representer, MonCode
from src.build_output.draw_hmm import draw_hmm
from src.generic.svg import svg_with_label, join_svgs_in_rectangle
import yaml
from pathlib import Path
from more_itertools import split_when

class HMM(NamedTuple):
    bgc_variant_id: BGC_Variant_ID
    transitions: List[List[Tuple[StateIdx, LogProb]]]  # u -> [(v, log_prob(u -> v))]
    emissions: List[List[LogProb]]  # u -> [log_prob(u -> emission)]

    state_types: List[DetailedHMMStateType]
    module_match_states: List[StateIdx]
    module_start_states: List[StateIdx]
    module_names: List[str]

    def to_json_dict(self):
        _emissions = [[emission if emission != LogProb(float('-inf')) else None
                      for emission in state_emissions]
                     for state_emissions in self.emissions]
        _transitions = [
            [
                (v, (lp if lp != LogProb(float('-inf')) else None))
                for v, lp in state_transitions
            ]
            for state_transitions in self.transitions
        ]
        return {
            'bgc_variant_id': self.bgc_variant_id.to_dict(),
            'transitions': _transitions,
            'emissions': _emissions,
            'state_types': [
                self.state_types[state_idx].name
                for state_idx in range(len(self.state_types))
            ],
            'module_match_states': self.module_match_states,
            'module_start_states': self.module_start_states,
            'module_names': self.module_names,
        }

    def from_json_dict(json_data: dict) -> HMM:
        return HMM(
            bgc_variant_id=BGC_Variant_ID.from_dict(json_data['bgc_variant_id']),
            transitions=[
                [
                    (v, (lp if lp is not None else LogProb(float('-inf'))))
                    for v, lp in state_transitions
                ]
                for state_transitions in json_data['transitions']
            ],
            emissions=[
                [emission if emission is not None else LogProb(float('-inf'))
                 for emission in state_emissions]
                for state_emissions in json_data['emissions']
            ],
            state_types=[
                DetailedHMMStateType[state_type]
                for state_type in json_data['state_types']
            ],
            module_start_states=json_data['module_start_states'],
            module_match_states=json_data['module_match_states'],
            module_names=json_data['module_names']
        )

    def is_emitting_state(self, state_idx: StateIdx) -> bool:
        return len(self.emissions[state_idx]) > 0

    def start_state_idx(self) -> StateIdx:
        return 0

    def draw(self,
             filename: Path,
             highlight_path: Optional[List[int]] = None,
             emission_names: Optional[List[str]] = None):
        if highlight_path is None:
            draw_hmm(self, filename)
            return

        if emission_names is not None:
            emission_labels = emission_names
        else:
            num_emitting_states = sum(1 for state_idx in highlight_path
                                      if self.is_emitting_state(state_idx))
            emission_labels = ['' for _ in range(num_emitting_states)]

        labels_iter = iter(emission_labels)
        whole_path_with_emissions = [
            (state_idx, next(labels_iter) if self.is_emitting_state(state_idx) else None)
            for state_idx in highlight_path
        ]

        # split the path into bgc iterations to avoid clutter
        # split *before* each new INITIAL state (except the very first element)
        paths_with_emissions = list(
            split_when(whole_path_with_emissions,
                        lambda _prev, curr: curr[0] == self.start_state_idx())
        )

        iteration_svgs = []
        for path_with_emissions in paths_with_emissions:
            path_part  = [state_idx for state_idx, _ in path_with_emissions]
            _ = draw_hmm(self, filename, path_part)
            figure_label = ', '.join(  
                (f'{state_idx}*{emission}' if emission is not None
                 else f'{state_idx}')
                for state_idx, emission in path_with_emissions
            )
            svg_text = filename.read_text()
            iteration_svgs.append(
                svg_with_label(
                    svg_text,
                    figure_label,
                    position='bottom',
                    font_size=18
                )
            )

        # combine iteration svgs into one svg
        svg_joined = join_svgs_in_rectangle([[iter_svg] for iter_svg in iteration_svgs])
            
        if emission_names is not None:
            title = '\n' + ', '.join(f'{emisson_name}' for emisson_name in emission_names)
            svg_joined = svg_with_label(svg_joined, title, position='top')

        filename.write_text(svg_joined)
        
        
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

