from typing import List, Tuple, NamedTuple, NewType
from collections import defaultdict
import math

StateIdx = int
LogProb = float
DiscreteProb = int

# DiscreteProb is a value between MIN_DISCRETE_VALUE and MAX_DISCRETE_PROB
MAX_DISCRETE_PROB = DiscreteProb(1000)
MIN_DISCRETE_PROB = DiscreteProb(0)

# LogProbs lesser than MIN_VALUE discretize to MIN_DISCRETE_PROB
THRESHOLD_MIN_VALUE = LogProb(-3)
# MIN_DISCRETE_PROB integrates to -inf
REDUNDANT_VALUE = float('-inf')


class HMM(NamedTuple):
    transitions: List[List[Tuple[StateIdx, LogProb]]]  # u -> [(v, log_prob(u -> v))]
    emissions: List[List[LogProb]]  # u -> [log_prob(u -> emission)]


def to_discrete_prob(prob: LogProb) -> DiscreteProb:
    if prob < THRESHOLD_MIN_VALUE:
        return MIN_DISCRETE_PROB
    else:
        # Map [THRESHOLD_MIN_VALUE, 0] to [MIN_DISCRETE_PROB, MAX_DISCRETE_PROB]
        disc_prob = round((prob - THRESHOLD_MIN_VALUE) /
                             (0 - THRESHOLD_MIN_VALUE) * (MAX_DISCRETE_PROB - MIN_DISCRETE_PROB) + MIN_DISCRETE_PROB)

        return DiscreteProb(disc_prob)


def to_prob(disc_prob: DiscreteProb) -> LogProb:
    if disc_prob == MIN_DISCRETE_PROB:
        return REDUNDANT_VALUE
    else:
        # Map [MIN_DISCRETE_PROB, MAX_DISCRETE_PROB] to [THRESHOLD_MIN_VALUE, 0]
        prob = (THRESHOLD_MIN_VALUE + (disc_prob - MIN_DISCRETE_PROB) /
                       (MAX_DISCRETE_PROB - MIN_DISCRETE_PROB) * (0 - THRESHOLD_MIN_VALUE))
        return LogProb(prob)


# TODO: apply for the list of all possible thresholds
def compute_p_value(hmm: HMM, threshold_score: LogProb) -> float:
    """
    computes p-value of given HMM with respect to the given threshold score
    as a sum of probabilities of all paths that have a score greater than the threshold
    multiplied by its quantity
    """
    disc_threshold_prob = to_discrete_prob(threshold_score)

    num_states = len(hmm.transitions)
    initial_state = 0
    final_state = num_states - 1

    # dp[state][discrete_prob] = number of paths
    dp = [defaultdict(int) for _ in range(num_states)]

    dp[initial_state][1000] = 1

    # Iterate over probabilities first, from highest to lowest
    for disc_prob in range(MAX_DISCRETE_PROB, disc_threshold_prob - 1, -1):
        # Then iterate over states
        for state in range(num_states):
            # Skip if no paths with this probability for this state
            if disc_prob not in dp[state]:
                continue

            number_of_paths = dp[state][disc_prob]
            path_prob = to_prob(disc_prob)

            # Process each outgoing transition
            for next_state, transition_prob in hmm.transitions[state]:
                # Skip if the transition probability is below the threshold
                if transition_prob < threshold_score:
                    continue

                # Process each emission of the next state if they exist
                if hmm.emissions[next_state]:
                    for emission_prob in hmm.emissions[next_state]:
                        # Skip if the emission probability is below the threshold
                        if emission_prob < threshold_score:
                            continue

                        new_path_prob = path_prob + transition_prob + emission_prob
                        disc_new_path_prob = to_discrete_prob(new_path_prob)

                        # Skip if the new probability is below the threshold
                        if disc_new_path_prob < disc_threshold_prob:
                            continue

                        dp[next_state][disc_new_path_prob] += number_of_paths
                else:
                    new_path_prob = path_prob + transition_prob
                    disc_new_path_prob = to_discrete_prob(new_path_prob)

                    # Skip if the new probability is below the threshold
                    if disc_new_path_prob < disc_threshold_prob:
                        continue

                    dp[next_state][disc_new_path_prob] += number_of_paths

    p_value = 0
    for disc_prob, number_of_paths in dp[final_state].items():
        if disc_prob >= disc_threshold_prob:
            # Yield p-value in range [0, 1]
            actual_prob = math.exp(to_prob(disc_prob))
            p_value += actual_prob * number_of_paths

    return p_value