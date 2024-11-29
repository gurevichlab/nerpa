from typing import List, Tuple, NamedTuple
from src.data_types import LogProb
#from src.matching.matcher_viterbi_types import HMM


class HMM(NamedTuple):
    adj_list: List[List[Tuple[int, LogProb]]]  # u -> [(v, log_prob(u -> v))]
    emission_log_probs: List[List[LogProb]]  # u -> [log_prob(u -> emission)]


def viterbi_algorithm(hmm: HMM, seq: List[int]) -> List[int]:
    """
    Viterbi algorithm for HMMs with detailed states and edges.

    Args:
    hmm: HMM - the HMM to run the Viterbi algorithm on
    seq: List[int] - the sequence of emissions

    Returns:
    List[int] - the most likely path through the HMM
    """
    n = len(seq)  # length of the sequence
    m = len(hmm.adj_list)  # number of states
    initial_state = 0
    final_state = m - 1

    # TODO: test if queue is faster than basic dynamic programming
    q = [(initial_state, 0, 0)]  # (state, num symbols, log probability)
    dp = [[float('-inf')] * (n + 1) for _ in range(m)]  # dp[u][i] = max log probability of a path of ending in state u with i symbols processed
    prev = [[-1] * (n + 1) for _ in range(m)]  # prev[u][i] = previous state of state u if i symbols have been processed
    while q:
        u, i, log_prob = q.pop()
        if i == n and hmm.adj_list[u]:  # if u should emit a symbol but we have processed all symbols
            continue
        for v, edge_log_prob in hmm.adj_list[u]:
            if hmm.emission_log_probs[v]:
                new_log_prob = log_prob + edge_log_prob + hmm.emission_log_probs[v][seq[i]]
                num_symbols = i + 1
            else:
                new_log_prob = log_prob + edge_log_prob
                num_symbols = i
            if new_log_prob > dp[v][num_symbols]:
                q.append((v, num_symbols, new_log_prob))
                prev[num_symbols][v] = u


    # Reconstruct the most likely path
    ext_path = [(final_state, n)]
    u, k = final_state, n
    while u != initial_state:
        u, k = prev[u][k], k - int(bool(hmm.emission_log_probs[u]))
        ext_path.append((u, k))
    return [state for state, _ in reversed(ext_path)]