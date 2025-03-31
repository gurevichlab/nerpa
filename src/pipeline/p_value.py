from typing import List, Tuple, NamedTuple, NewType
import math
import numpy as np

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
    def __init__(self, hmm: HMM):
        self.hmm = hmm
        self.p_values = compute_hmm_p_values(hmm)

    def __call__(self, threshold_score: LogProb) -> Prob:
        disc_threshold_score = to_discrete_prob(threshold_score)
        p_value = self.p_values[disc_threshold_score]
        return p_value


def to_discrete_prob(prob: LogProb) -> DiscreteProb:
    if prob < LOG_PROB_MIN_VALUE:
        return MIN_DISCRETE_PROB
    else:
        # Map [LOG_PROB_MIN_VALUE, LOG_PROB_MAX_VALUE] to [MIN_DISCRETE_PROB, MAX_DISCRETE_PROB]
        disc_prob = round(MIN_DISCRETE_PROB + (prob - LOG_PROB_MIN_VALUE) /
                            (LOG_PROB_MAX_VALUE - LOG_PROB_MIN_VALUE) * (MAX_DISCRETE_PROB - MIN_DISCRETE_PROB))

        return DiscreteProb(disc_prob)


def to_prob(disc_prob: DiscreteProb) -> LogProb:
    if disc_prob == MIN_DISCRETE_PROB:
        return float('-inf')
    else:
        # Map [MIN_DISCRETE_PROB, MAX_DISCRETE_PROB] to [LOG_PROB_MIN_VALUE, LOG_PROB_MAX_VALUE]
        log_prob = (LOG_PROB_MIN_VALUE + (disc_prob - MIN_DISCRETE_PROB) /
                       (MAX_DISCRETE_PROB - MIN_DISCRETE_PROB) * (LOG_PROB_MAX_VALUE - LOG_PROB_MIN_VALUE))
        return LogProb(log_prob)


def compute_hmm_p_values(hmm: HMM) -> np.ndarray:
    num_states = len(hmm.transitions)
    initial_state = 0
    final_state = num_states - 1

    # dp[state, discrete_prob] = number of paths ending at the state 'state' and
    # having discretized probability 'discrete_prob'
    dp = np.zeros((num_states, MAX_DISCRETE_PROB + 1), dtype=int)

    dp[initial_state][MAX_DISCRETE_PROB] = 1

    # Iterate over probabilities first, from highest to lowest
    for disc_prob in range(MAX_DISCRETE_PROB, MIN_DISCRETE_PROB - 1, -1):
        # Then iterate over states
        for state in range(num_states):
            # Skip if no paths with this probability for this state
            if dp[state, disc_prob] == 0:
                continue

            number_of_paths = dp[state, disc_prob]
            path_prob = to_prob(DiscreteProb(disc_prob))

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

    p_values = np.zeros((MAX_DISCRETE_PROB - MIN_DISCRETE_PROB + 1), dtype=float)

    for k in range(MIN_DISCRETE_PROB, MAX_DISCRETE_PROB + 1):
        p_value = 0.0
        for disc_prob in range(k, MAX_DISCRETE_PROB + 1):
            actual_prob = math.exp(to_prob(DiscreteProb(disc_prob)))
            p_value += actual_prob * dp[final_state, disc_prob]

        # Not count the contribution of the found path itself
        p_value -= math.exp(to_prob(DiscreteProb(k)))

        # Prevent possible miscalculation caused by discretization
        p_value = min(max(p_value, 0.0), 1.0)

        p_values[k] = Prob(p_value)

    return p_values