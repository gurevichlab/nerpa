from __future__ import annotations

from src.hmm.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    DetailedHMMStateType,
    DetailedHMMState,
    StateIdx,
)

ET = DetailedHMMEdgeType
ST = DetailedHMMStateType


def add_skips_at_end(constructor):  # constructor: HMM_Constructor
    def is_last_module_state(state_idx: StateIdx) -> bool:
        return constructor.states[state_idx].related_module_idx == len(constructor.bgc_variant.modules) - 1

    matching_finished_states = [
        state_idx
        for state_idx, state in enumerate(constructor.states)
        if (state.state_type == ST.MATCHING_FINISHED and
            not is_last_module_state(state_idx))  # no skips after the last module
    ]
    matching_finished_states = sorted(matching_finished_states,
                                      key=lambda idx: constructor.states[idx].related_module_idx)


    next_skipping_at_end_state_idx = constructor.final_state
    for matching_finished_state in reversed(matching_finished_states):
        module_idx = constructor.states[matching_finished_state].related_module_idx
        module = constructor.bgc_variant.modules[module_idx]

        constructor.states.extend([
            DetailedHMMState(state_type=ST.INSERT_AT_END,
                             emissions=constructor.get_insert_emissions(module),
                             related_module_idx=module_idx),
            DetailedHMMState(state_type=ST.SKIPPING_MODULES_AT_END,
                             related_module_idx=module_idx + 1)
        ])

        insert_at_end = len(constructor.states) - 2
        skipping_at_end = len(constructor.states) - 1

        constructor.add_edge(matching_finished_state, insert_at_end, ET.START_INSERTING)
        constructor.add_edge(matching_finished_state, skipping_at_end , ET.NO_INSERTIONS)

        constructor.add_edge(insert_at_end, insert_at_end, ET.INSERT_AT_END)
        constructor.add_edge(insert_at_end, skipping_at_end, ET.FINISH_INSERTING)

        constructor.add_edge(skipping_at_end, next_skipping_at_end_state_idx, ET.SKIP_MODULES_AT_END)
        next_skipping_at_end_state_idx = skipping_at_end







