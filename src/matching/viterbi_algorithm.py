from typing import List, Tuple, Optional
from src.data_types import LogProb
from src.matching.hmm_auxiliary_types import HMM, StateIdx
from src.monomer_names_helper import MonCode
from src.data_types import LogProb
import numpy as np
from math import log


# checkpoint: (state, symbol_count) means that
# at state `state` the symbol observed_sequence[symbol_count] was emitted
def get_opt_path_with_score(hmm: HMM,
                            observed_sequence: List[MonCode],
                            checkpoints: Optional[List[Tuple[StateIdx, int]]] = None) -> Tuple[LogProb, List[StateIdx]]:
    seq_len = len(observed_sequence)
    num_states = len(hmm.transitions)

    if checkpoints is None:
        checkpoints = [(0, 0), (num_states - 1, seq_len)]

    # Initialize DP table for log probabilities
    dp = np.full((num_states, seq_len + 1), -np.inf)
    prev = np.full((num_states, seq_len + 1), -1, dtype=int)

    # Start from the first checkpoint
    first_state, first_symbols = checkpoints[0]
    dp[first_state][first_symbols] = 0.0

    # Forward pass: Process all checkpoints
    for chk_idx in range(len(checkpoints) - 1):
        cur_state, cur_sym_count = checkpoints[chk_idx]
        next_state, next_sym_count = checkpoints[chk_idx + 1]

        for sym_count in range(cur_sym_count, next_sym_count + 1):
            for state in range(cur_state, next_state):
                if sym_count == next_sym_count and hmm.emissions[state]:
                    continue

                log_prob = dp[state][sym_count]
                if log_prob == -np.inf:
                    continue

                for edge_to, edge_log_prob in hmm.transitions[state]:
                    new_sym_count = sym_count
                    new_log_prob = log_prob + edge_log_prob

                    if hmm.emissions[state]:  # If emitting, add emission probability
                        new_log_prob += hmm.emissions[state][observed_sequence[sym_count]]
                        new_sym_count += 1

                    if new_log_prob > dp[edge_to][new_sym_count]:
                        dp[edge_to][new_sym_count] = new_log_prob
                        prev[edge_to][new_sym_count] = state

    # Extract best final probability
    final_state, final_symbols = checkpoints[-1]
    best_score = dp[final_state][final_symbols]

    # Backtrace to reconstruct the optimal path
    path = [final_state]
    state = final_state
    symbols = final_symbols

    while state != first_state:
        prev_state = prev[state][symbols]
        if prev_state == -1:
            raise RuntimeError("Backtrace failed: No valid path found!")

        if hmm.emissions[prev_state]:  # If emitting, move symbol index
            symbols -= 1

        state = prev_state
        path.append(state)

    path.reverse()

    return float(best_score), path

