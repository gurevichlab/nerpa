from typing import List, Tuple, NamedTuple, NewType
import numpy as np
from numpy.typing import NDArray

StateIdx = NewType('StateIdx', int)
LogProb = NewType('LogProb', float | float('-inf'))
DiscreteProb = NewType('DiscreteProb', int)
Prob = NewType('Prob', float)  # actual probability value in range [0, 1]

# DiscreteProb is a value between MIN_DISCRETE_VALUE and MAX_DISCRETE_PROB
MAX_DISCRETE_PROB = DiscreteProb(1000)
MIN_DISCRETE_PROB = DiscreteProb(0)

# LogProbs less than MIN_VALUE discretize to MIN_DISCRETE_PROB
LOG_PROB_MIN_VALUE = LogProb(-30)
LOG_PROB_MAX_VALUE = LogProb(0)


class HMM(NamedTuple):
    transitions: List[List[Tuple[StateIdx, LogProb]]]   # u -> [(v, log_prob(u -> v))]
    emissions: List[List[LogProb]]  # u -> [log_prob(u -> emission)]


class PValueEstimator:
    _p_values: NDArray[np.float64]

    def __init__(self, hmm: HMM):
        self._p_values = compute_hmm_p_values(hmm)

    def __call__(self, path_log_prob: LogProb) -> Prob:
        disc_threshold_score = to_discrete_prob(path_log_prob)
        p_value = self._p_values[disc_threshold_score]
        return p_value


def to_discrete_prob(log_prob: LogProb) -> DiscreteProb:
    if log_prob < LOG_PROB_MIN_VALUE:
        return MIN_DISCRETE_PROB
    else:
        # Map [LOG_PROB_MIN_VALUE, LOG_PROB_MAX_VALUE] to [MIN_DISCRETE_PROB, MAX_DISCRETE_PROB]
        disc_prob = round(MIN_DISCRETE_PROB + (log_prob - LOG_PROB_MIN_VALUE) /
                            (LOG_PROB_MAX_VALUE - LOG_PROB_MIN_VALUE) * (MAX_DISCRETE_PROB - MIN_DISCRETE_PROB))

        return DiscreteProb(disc_prob)


def to_log_prob(disc_prob: DiscreteProb) -> LogProb:
    if disc_prob == MIN_DISCRETE_PROB:
        return LogProb(float('-inf'))
    else:
        # Map [MIN_DISCRETE_PROB, MAX_DISCRETE_PROB] to [LOG_PROB_MIN_VALUE, LOG_PROB_MAX_VALUE]
        log_prob = (LOG_PROB_MIN_VALUE + (disc_prob - MIN_DISCRETE_PROB) /
                       (MAX_DISCRETE_PROB - MIN_DISCRETE_PROB) * (LOG_PROB_MAX_VALUE - LOG_PROB_MIN_VALUE))
        return LogProb(log_prob)


def compute_hmm_p_values(hmm: HMM) -> NDArray[np.float64]:
    """
    Returns an array p_values of length (MAX_DISCRETE_PROB - MIN_DISCRETE_PROB + 1).
    For each k, the element p_values[k] is the probability of finding a path with
    the score to_log_prob(k) or higher in the HMM. One path with exactly that score
    is not counted in the probability.
    """
    num_states = len(hmm.transitions)
    initial_state = 0
    final_state = num_states - 1

    # dp[state, discrete_prob] = number of paths ending at the state 'state' and
    # having discretized probability 'discrete_prob'
    dp = np.zeros((num_states, MAX_DISCRETE_PROB + 1), dtype=int)

    dp[initial_state, MAX_DISCRETE_PROB] = 1

    # Iterate over probabilities first, from highest to lowest
    for disc_prob in map(DiscreteProb, range(MAX_DISCRETE_PROB, MIN_DISCRETE_PROB - 1, -1)):
        # Then iterate over states
        for state in range(num_states):
            # Skip if no paths with this probability for this state
            if dp[state, disc_prob] == 0:
                continue

            number_of_paths = dp[state, disc_prob]
            path_prob = to_log_prob(disc_prob)

            # Process each outgoing transition
            for next_state, transition_prob in hmm.transitions[state]:
                # Process each emission of the next state if they exist
                if hmm.emissions[next_state]:
                    for emission_prob in hmm.emissions[next_state]:
                        new_path_prob = path_prob + transition_prob + emission_prob
                        disc_new_path_prob = to_discrete_prob(new_path_prob)

                        dp[next_state][disc_new_path_prob] += number_of_paths
                else:
                    new_path_prob = path_prob + transition_prob
                    disc_new_path_prob = to_discrete_prob(new_path_prob)

                    dp[next_state][disc_new_path_prob] += number_of_paths

    def to_prob(disc_prob: DiscreteProb) -> np.float64:
        return np.exp(to_log_prob(DiscreteProb(disc_prob)))

    p_values = np.zeros((MAX_DISCRETE_PROB - MIN_DISCRETE_PROB + 1), dtype=np.float64)

    p_values[MAX_DISCRETE_PROB] = to_prob(MAX_DISCRETE_PROB) * dp[final_state, MAX_DISCRETE_PROB]

    # Iterate over probabilities from highest to lowest
    for disc_prob in map(DiscreteProb, range(MAX_DISCRETE_PROB - 1, MIN_DISCRETE_PROB - 1, -1)):
        p_values[disc_prob] = p_values[disc_prob + 1] + to_prob(disc_prob) * dp[final_state, disc_prob]

    # Not count the contribution of the found path itself
    contributions = np.array([to_prob(k) for k in range(MIN_DISCRETE_PROB, MAX_DISCRETE_PROB + 1)])
    p_values -= contributions

    # Prevent possible miscalculation caused by discretization
    p_values = np.clip(p_values, 0.0, 1.0)

    return p_values