from __future__ import annotations

import math
from math import e
from typing import List, Tuple, Dict, NamedTuple, NewType
from src.data_types import LogProb, Prob
from src.matching.hmm_auxiliary_types import HMM, HMM_LOUC, StateIdx, HMM_LPKC
import numpy as np
from numpy.typing import NDArray
from joblib import Parallel, delayed

DiscreteLogOdds = NewType('DiscreteLogOdds', int)

# DiscreteProb is a value between MIN_DISCRETE_VALUE and MAX_DISCRETE_PROB
MIN_DISCRETE_LOG_ODDS = DiscreteLogOdds(0)
MAX_DISCRETE_LOG_ODDS = DiscreteLogOdds(50000)

# LogProbs less than MIN_VALUE discretize to MIN_DISCRETE_PROB
LOG_ODDS_MIN_VALUE = -20.0
LOG_ODDS_MAX_VALUE = 50.0


class PValueEstimator:
    _p_values: NDArray[np.float64]

    def __init__(self,
                 hmm_lp: HMM_LPKC,
                 hmm_lo: HMM_LOUC):
        self._p_values = compute_hmm_p_values(hmm_lp, hmm_lo)

    def __call__(self,
                 path_log_odds: float,
                 path_log_prob: LogProb) -> Prob:
        disc_threshold_score = to_discrete_log_odds(path_log_odds)
        p_value = float(self._p_values[disc_threshold_score]) - math.e ** path_log_prob
        return Prob(max(p_value, 0))

    @classmethod
    def _from_precomputed_p_values(cls, p_values: NDArray[np.float64]) -> PValueEstimator:
        """
        Create a PValueEstimator from precomputed p-values.
        """
        obj = cls.__new__(cls)
        obj._p_values = p_values
        return obj

    @classmethod
    def _precompute_p_values_for_hmms(cls,
                                      hmms: List[Tuple[HMM_LPKC, HMM_LOUC]],
                                      num_threads: int = 1) -> List[NDArray[np.float64]]:
        p_values_list = Parallel(n_jobs=num_threads, backend="loky")(
            delayed(compute_hmm_p_values)(hmm_lp, hmm_lo) for hmm_lp, hmm_lo in hmms
        )
        return p_values_list

    @classmethod
    def precompute_p_value_estimators_for_hmms(cls,
                                               hmms: List[Tuple[HMM_LPKC, HMM_LOUC]],
                                               num_threads: int = 1):
        """
        Precompute p-values for a list of HMMs in parallel and return a list of PValueEstimator objects.
        This method is used to speed up p-value estimation for multiple HMMs.
        """
        p_values_list = cls._precompute_p_values_for_hmms(hmms, num_threads)
        return [PValueEstimator._from_precomputed_p_values(p_values) for p_values in p_values_list]


def to_discrete_log_odds(log_odds: float) -> DiscreteLogOdds:
    if log_odds < LOG_ODDS_MIN_VALUE:
        return MIN_DISCRETE_LOG_ODDS
    elif log_odds > LOG_ODDS_MAX_VALUE:
        return MAX_DISCRETE_LOG_ODDS
    else:
        # Map [LOG_PROB_MIN_VALUE, LOG_PROB_MAX_VALUE] to [MIN_DISCRETE_PROB, MAX_DISCRETE_PROB]
        disc_log_odds = round(MIN_DISCRETE_LOG_ODDS + (log_odds - LOG_ODDS_MIN_VALUE) /
                          (LOG_ODDS_MAX_VALUE - LOG_ODDS_MIN_VALUE) * (MAX_DISCRETE_LOG_ODDS - MIN_DISCRETE_LOG_ODDS))

        return DiscreteLogOdds(disc_log_odds)


def to_log_odds(disc_log_odds: DiscreteLogOdds) -> float:
    if disc_log_odds == MIN_DISCRETE_LOG_ODDS:
        return float('-inf')
    else:
        # Map [MIN_DISCRETE_PROB, MAX_DISCRETE_PROB] to [LOG_PROB_MIN_VALUE, LOG_PROB_MAX_VALUE]
        log_odds = (LOG_ODDS_MIN_VALUE + (disc_log_odds - MIN_DISCRETE_LOG_ODDS) /
                    (MAX_DISCRETE_LOG_ODDS - MIN_DISCRETE_LOG_ODDS) * (LOG_ODDS_MAX_VALUE - LOG_ODDS_MIN_VALUE))
        return log_odds


def compute_hmm_p_values(hmm_lp: HMM_LPKC,
                         hmm_lo: HMM_LOUC) -> NDArray[np.float64]:
    """
    Returns an array p_values of length (MAX_DISCRETE_PROB - MIN_DISCRETE_PROB + 1).
    For each k, the element p_values[k] is the probability of finding a path with
    the score to_log_prob(k) or higher in the HMM. One path with exactly that score
    is not counted in the probability.
    """
    eps = 1e-5  # Small value to avoid numerical issues
    u = next((u
             for u, edges in enumerate(hmm_lp.transitions)
             if edges and abs(sum(e ** transition_lp for _, transition_lp in edges) - 1) > eps),
            None)
    if u is not None:
        raise ValueError(f"HMM {hmm_lp.bgc_variant_id} invalid: transition probabilities for the state {u} do not sum to 1.")
    u = next((u
              for u, emissions in enumerate(hmm_lp.emissions)
              if emissions and abs(sum(e ** emission_lp for emission_lp in emissions) - 1) > eps),
             None)
    if u is not None:
        raise ValueError(f"HMM {hmm_lp.bgc_variant_id} invalid: emission probabilities for the state {u} do not sum to "
                         f"{sum(e ** emission_lp for emission_lp in hmm_lp.emissions[u])} instead of 1.")

    num_states = len(hmm_lp.transitions)
    initial_state = 0
    final_state = num_states - 1

    # dp[state, discrete_prob] = total probability of paths ending at the state 'state' and
    # having discretized log-odds score 'discrete_prob'
    dp = np.zeros((num_states, MAX_DISCRETE_LOG_ODDS + 1), dtype=np.float64)

    dp[initial_state, to_discrete_log_odds(0.0)] = 1.0

    # Iterate over probabilities first, from highest to lowest
    # MIN_DISCRETE_PROB is excluded to avoid overflow: it corresponds to -inf log probability and infinite number of paths
    for state in range(num_states):
        for disc_log_odds in map(DiscreteLogOdds, range(MAX_DISCRETE_LOG_ODDS, MIN_DISCRETE_LOG_ODDS, -1)):
            # Skip if no paths with this probability for this state
            if dp[state, disc_log_odds] == 0.0:
                continue

            assert 0.0 <= dp[state, disc_log_odds] <= 1.0, \
                f"Invalid probability {dp[state, disc_log_odds]} for state {state} and discrete probability {disc_log_odds}"

            current_prob = dp[state, disc_log_odds]
            current_lo = to_log_odds(disc_log_odds)
            # Process each outgoing transition
            for i, (next_state, transition_lp) in enumerate(hmm_lp.transitions[state]):
                next_state_lo, transition_lo = hmm_lo.transitions[state][i]
                assert next_state_lo == next_state, \
                    f"Transition from state {state} to {next_state} has inconsistent log-odds and log-probabilities: " \

                # Process each emission of the next state if they exist
                if hmm_lp.emissions[next_state]:
                    for mon_code, emission_lp in enumerate(hmm_lp.emissions[next_state]):
                        emission_lo = hmm_lo.emissions[next_state][mon_code]
                        new_lo = current_lo + transition_lo + emission_lo
                        new_prob = current_prob * math.e ** transition_lp * math.e ** emission_lp
                        disc_new_lo = to_discrete_log_odds(new_lo)

                        dp[next_state][disc_new_lo] += new_prob
                else:
                    new_lo = current_lo + transition_lo
                    new_prob = current_prob * math.e ** transition_lp
                    disc_new_lo = to_discrete_log_odds(new_lo)

                    dp[next_state][disc_new_lo] += new_prob

    p_values = np.zeros((MAX_DISCRETE_LOG_ODDS - MIN_DISCRETE_LOG_ODDS + 1), dtype=np.float64)

    p_values[MAX_DISCRETE_LOG_ODDS] = dp[final_state, MAX_DISCRETE_LOG_ODDS]

    # Iterate over probabilities from highest to lowest
    for disc_log_odds in map(DiscreteLogOdds, range(MAX_DISCRETE_LOG_ODDS - 1, MIN_DISCRETE_LOG_ODDS - 1, -1)):
        p_values[disc_log_odds] = p_values[disc_log_odds + 1] + dp[final_state, disc_log_odds]

    # Prevent possible miscalculation caused by discretization
    p_values = np.clip(p_values, 0.0, 1.0)

    return p_values