from typing import List, Tuple, NamedTuple
from src.data_types import LogProb
from src.matching.hmm_auxiliary_types import HMM
from src.data_types import LogProb
from math import log


def get_opt_path_with_emissions(hmm: HMM,
                                initial_state: int,
                                final_state: int,
                                observed_sequence: List[int]) -> Tuple[LogProb, List[Tuple[int, int]]]:
    n = len(observed_sequence)  # length of the sequence
    m = len(hmm.adj_list)  # number of states

    q = [(initial_state, 0, 0)]  # (state, num symbols, log probability)
    dp = [[float('-inf')] * (n + 1)
          for _ in range(m)]  # dp[u][i] = max log probability for a path of ending in state u with i symbols processed
    prev = [[-1] * (n + 1) for _ in range(m)]  # prev[u][i] = previous state of state u if i symbols have been processed

    while q:
        u, i, log_prob = q.pop(0)
        if i == n and hmm.emission_log_probs[u]:  # if u should emit a symbol but we have processed all symbols
            continue
        for v, edge_log_prob in hmm.adj_list[u]:
            if hmm.emission_log_probs[u]:  # if u emits a symbol
                new_log_prob = log_prob + edge_log_prob + hmm.emission_log_probs[u][observed_sequence[i]]
                num_symbols = i + 1
            else:
                new_log_prob = log_prob + edge_log_prob
                num_symbols = i
            if new_log_prob > dp[v][num_symbols]:
                dp[v][num_symbols] = new_log_prob
                prev[v][num_symbols] = u
                q.append((v, num_symbols, new_log_prob))

    # Reconstruct the most likely path
    # final_state = max([i for i in range(m) if emission_matrix[i][observed_sequence[n-1]] > 0],
    #                  key=lambda k: dp[k][n-1] + log(emission_matrix[k][observed_sequence[n-1]]))
    ext_path = [(final_state, n)]
    emitted_symbols = [(final_state, None)]
    u, k = final_state, n
    while u != initial_state:
        prev_state = prev[u][k]
        assert prev_state != -1, 'prev_state not defined!!!'
        symbol_was_emitted = bool(hmm.emission_log_probs[prev_state])
        emitted_symbol = observed_sequence[k - 1] if symbol_was_emitted else None
        u, k = prev_state, k - int(symbol_was_emitted)
        ext_path.append((u, k))
        emitted_symbols.append((u, emitted_symbol))

    path = [state for state, symbol in reversed(emitted_symbols)]
    assert path[0] == initial_state and path[-1] == final_state, 'the path endpoints are wrong!!!'
    return dp[final_state][n], list(reversed(emitted_symbols))


def get_hmm_score(hmm: HMM,
                  observed_sequence: List[int]) -> LogProb:
    n = len(observed_sequence)  # length of the sequence
    m = len(hmm.adj_list)  # number of states

    initial_state = 0
    final_state = m - 1
    q = [(initial_state, 0, 0)]  # (state, num symbols, log probability)
    dp = [[float('-inf')] * (n + 1)
          for _ in range(m)]  # dp[u][i] = max log probability for a path of ending in state u with i symbols processed

    while q:
        u, i, log_prob = q.pop(0)
        if i == n and hmm.emission_log_probs[u]:  # if u should emit a symbol but we have processed all symbols
            continue
        for v, edge_log_prob in hmm.adj_list[u]:
            if hmm.emission_log_probs[u]:  # if u emits a symbol
                new_log_prob = log_prob + edge_log_prob + hmm.emission_log_probs[u][observed_sequence[i]]
                num_symbols = i + 1
            else:
                new_log_prob = log_prob + edge_log_prob
                num_symbols = i
            if new_log_prob > dp[v][num_symbols]:
                dp[v][num_symbols] = new_log_prob
                q.append((v, num_symbols, new_log_prob))

    return dp[final_state][n]
