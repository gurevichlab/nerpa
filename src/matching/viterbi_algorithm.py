from typing import List, Tuple, Optional
from src.data_types import LogProb
from src.matching.hmm_auxiliary_types import HMM, StateIdx
from src.monomer_names_helper import MonCode
from src.data_types import LogProb
import numpy as np
from math import log


def get_opt_path_with_score(hmm: HMM,
                            observed_sequence: List[MonCode],
                            start_state: Optional[StateIdx] = None,
                            finish_state: Optional[StateIdx] = None) -> Tuple[LogProb, List[StateIdx]]:
    seq_len = len(observed_sequence)
    num_states = len(hmm.transitions)

    if start_state is None:
        start_state = 0
    if finish_state is None:
        finish_state = num_states - 1

    # Initialize DP table for log probabilities
    dp = np.full((num_states, seq_len + 1), -np.inf)
    prev = np.full((num_states, seq_len + 1), -1, dtype=int)

    # Forward pass: Process all checkpoints
    for sym_count in range(len(observed_sequence) + 1):
        '''
        If we ignore gene/module iteration edges, states are sorted topologically.
        After an iteration edge, however, the state index is less than the previous one.
        Between any two iterations, a symbol should be emitted (otherwise iteration is useless).
        So I iterate over all states in the range [cur_state, next_state) and if an iteration edge is applied, 
        I iterate one more time over the same range. This is done at most once.
        '''
        iteration_applied = False
        for _dp_run in range(2):  # possibly go through states once more if iteration was applied
            if _dp_run == 1 and not iteration_applied:
                break

            for state in range(num_states):  # iterate over all states (not just [start_state, finish_state) because of iteration)
                if sym_count == len(observed_sequence) and hmm.emissions[state]:
                    continue

                log_prob = dp[state][sym_count]
                if (state, sym_count) != (start_state, 0) and prev[state][sym_count] == -1:  # unreachable state
                    continue

                for edge_to, edge_log_prob in hmm.transitions[state]:
                    new_sym_count = sym_count
                    new_log_prob = log_prob + edge_log_prob

                    if hmm.emissions[state]:  # If emitting, add emission probability
                        new_log_prob += hmm.emissions[state][observed_sequence[sym_count]]
                        new_sym_count += 1

                    if new_log_prob > dp[edge_to][new_sym_count] or prev[edge_to][new_sym_count] == -1:
                        dp[edge_to][new_sym_count] = new_log_prob
                        prev[edge_to][new_sym_count] = state
                        if edge_to < state and new_sym_count == sym_count:  # gene/module iteration
                            iteration_applied = True

    # Extract best final probability
    best_score = dp[finish_state][len(observed_sequence)]

    # Backtrace to reconstruct the optimal path
    path = [finish_state]
    state = finish_state
    symbols = len(observed_sequence)

    while state != start_state:
        prev_state = prev[state][symbols]
        if prev_state == -1:
            raise RuntimeError("Backtrace failed: No valid path found!")

        if hmm.emissions[prev_state]:  # If emitting, move symbol index
            symbols -= 1

        state = prev_state
        path.append(state)

    path.reverse()
    return float(best_score), path
