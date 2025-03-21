from typing import List, Tuple, NamedTuple
from collections import defaultdict
import math

StateIdx = int
LogProb = float


class HMM(NamedTuple):
    transitions: List[List[Tuple[StateIdx, LogProb]]]  # u -> [(v, log_prob(u -> v))]
    emissions: List[List[LogProb]]  # u -> [log_prob(u -> emission)]


# TODO: apply for the list of all possible thresholds
def compute_p_value(hmm: HMM, threshold_score: LogProb) -> float:
    """
    computes p-value of given HMM with respect to the given threshold score
    as a sum of probabilities of all paths that have a score greater than the threshold
    multiplied by its quantity
    """
    # Convert threshold log probability to actual probability, then discretize
    threshold_prob = math.exp(threshold_score)
    disc_threshold_prob = round(threshold_prob * 1000)

    num_states = len(hmm.transitions)
    initial_state = 0
    final_state = num_states - 1

    # dp[state][discrete_prob] = number of paths
    dp = [defaultdict(int) for _ in range(num_states)]

    dp[initial_state][1000] = 1

    # Process states in order
    for state in range(num_states):
        # Skip empty states
        if not dp[state]:
            continue

        # Process each probability in this state
        for disc_path_prob, number_of_paths in list(dp[state].items()):
            # Skip if the current path probability is below the threshold
            if disc_path_prob < disc_threshold_prob:
                continue

            path_prob = disc_path_prob / 1000.0

            # Process each outgoing transition
            for next_state, log_transition_prob in hmm.transitions[state]:
                # Skip if the transition probability is below the threshold
                if log_transition_prob < threshold_score:
                    continue

                transition_prob = math.exp(log_transition_prob)

                # Process each emission of the next state if they exist
                if hmm.emissions[next_state]:
                    for log_emission_prob in hmm.emissions[next_state]:
                        # Skip if the emission probability is below the threshold
                        if log_emission_prob < threshold_score:
                            continue

                        emission_prob = math.exp(log_emission_prob)

                        new_path_prob = path_prob * transition_prob * emission_prob
                        disc_new_path_prob = round(new_path_prob * 1000)

                        dp[next_state][disc_new_path_prob] += number_of_paths
                else:
                    new_path_prob = path_prob * transition_prob
                    disc_new_path_prob = round(new_path_prob * 1000)

                    dp[next_state][disc_new_path_prob] += number_of_paths

            if state != final_state:
                dp[state].clear()

    p_value = 0
    for disc_prob, number_of_paths in dp[final_state].items():
        actual_prob = disc_prob / 1000
        if actual_prob >= threshold_prob:
            p_value += actual_prob * number_of_paths

    return p_value