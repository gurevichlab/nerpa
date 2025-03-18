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
    assert checkpoints is None

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
            '''
            If we ignore gene/module iteration edges, states are sorted topologically.
            After an iteration edge, however, the state index is less than the previous one.
            Between any two iterations, a symbol should be emitted (otherwise iteration is useless).
            So I iterate over all states in the range [cur_state, next_state) and if an iteration edge is applied, 
            I iterate one more time over the same range. This is done at most once.
            '''
            iteration_applied = False
            for _dp_run in range(2):  # possibly go through range(cur_state, next_state) once more if iteration was applied
                if _dp_run == 1 and not iteration_applied:
                    break

                for state in range(cur_state, next_state):
                    if sym_count == next_sym_count and hmm.emissions[state]:
                        continue

                    log_prob = dp[state][sym_count]
                    if (state, sym_count) != checkpoints[chk_idx] and prev[state][sym_count] == -1:
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
                            if edge_to < state and next_sym_count == sym_count:  # gene/module iteration
                                iteration_applied = True

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

    check_score = 0.0
    mon_idx = 0
    for i in range(1, len(path)):
        check_score += next(weight
                            for (u, weight) in hmm.transitions[path[i-1]]
                            if u == path[i])
        if hmm.emissions[path[i-1]]:
            check_score += hmm.emissions[path[i-1]][observed_sequence[mon_idx]]
            mon_idx += 1

    if abs(check_score - best_score) > 1e-6:
        raise RuntimeError(f"Score mismatch: {check_score} != {best_score}")

    print(f"Score: {best_score}")
    print(f"Path: {path}")
    return float(best_score), path

