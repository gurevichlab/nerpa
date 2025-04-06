from typing import List, Tuple, Dict, NamedTuple, NewType
import numpy as np
from numpy.typing import NDArray

StateIdx = NewType('StateIdx', int)
LogProb = NewType('LogProb', float)
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
        return Prob(p_value)


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
            path_prob = to_log_prob(DiscreteProb(disc_prob))

            # Process each outgoing transition
            for next_state, transition_prob in hmm.transitions[state]:
                # Process each emission of the next state if they exist
                if hmm.emissions[next_state]:
                    for emission_prob in hmm.emissions[next_state]:
                        new_path_prob = path_prob + transition_prob + emission_prob
                        disc_new_path_prob = to_discrete_prob(LogProb(new_path_prob))

                        dp[next_state][disc_new_path_prob] += number_of_paths
                else:
                    new_path_prob = path_prob + transition_prob
                    disc_new_path_prob = to_discrete_prob(LogProb(new_path_prob))

                    dp[next_state][disc_new_path_prob] += number_of_paths

    def to_prob(disc_prob: DiscreteProb) -> np.float64:
        return np.exp(to_log_prob(DiscreteProb(disc_prob)))

    p_values = np.zeros((MAX_DISCRETE_PROB - MIN_DISCRETE_PROB + 1), dtype=np.float64)

    p_values[MAX_DISCRETE_PROB] = to_prob(MAX_DISCRETE_PROB) * dp[final_state, MAX_DISCRETE_PROB]

    # Iterate over probabilities from highest to lowest
    for disc_prob in map(DiscreteProb, range(MAX_DISCRETE_PROB - 1, MIN_DISCRETE_PROB - 1, -1)):
        p_values[disc_prob] = p_values[disc_prob + 1] + to_prob(DiscreteProb(disc_prob)) * dp[final_state, disc_prob]

    # Subtract the contribution of one path with each discrete probability
    contributions = np.array([to_prob(DiscreteProb(k))
                              for k in range(MIN_DISCRETE_PROB, MAX_DISCRETE_PROB + 1)])
    p_values -= contributions

    # Prevent possible miscalculation caused by discretization
    p_values = np.clip(p_values, 0.0, 1.0)

    return p_values


def get_all_paths_log_probs(hmm: HMM) -> List[LogProb]:
    path_log_probs = []
    num_states = len(hmm.transitions)
    initial_state = 0
    final_state = num_states - 1

    def path_log_prob_generator(state: StateIdx, current_log_prob: LogProb):
        if state == final_state:
            path_log_probs.append(current_log_prob)
            return

        for next_state, transition_prob in hmm.transitions[state]:
            new_log_prob = current_log_prob + transition_prob

            # Process each emission of the next state if they exist
            if hmm.emissions[next_state]:
                for emission_prob in hmm.emissions[next_state]:
                    path_log_prob_generator(next_state, LogProb(new_log_prob + emission_prob))
            else:
                path_log_prob_generator(next_state, LogProb(new_log_prob))

    # Start the recursion
    path_log_prob_generator(StateIdx(initial_state), LogProb(0))

    return path_log_probs


def compute_naive_p_values(hmm: HMM) -> NDArray[np.float64]:
    path_log_probs = get_all_paths_log_probs(hmm)
    path_disc_probs = [to_discrete_prob(log_prob) for log_prob in path_log_probs]

    path_counts = np.zeros(MAX_DISCRETE_PROB - MIN_DISCRETE_PROB + 1, dtype=int)
    for disc_prob in path_disc_probs:
        path_counts[disc_prob - MIN_DISCRETE_PROB] += 1

    p_values = np.zeros(MAX_DISCRETE_PROB - MIN_DISCRETE_PROB + 1, dtype=np.float64)

    def to_prob(disc_prob: int) -> float:
        return np.exp(to_log_prob(DiscreteProb(disc_prob)))

    p_values[MAX_DISCRETE_PROB - MIN_DISCRETE_PROB] = (to_prob(MAX_DISCRETE_PROB)
                                                       * path_counts[MAX_DISCRETE_PROB - MIN_DISCRETE_PROB])

    # Calculate p-values for the rest of the discrete probabilities
    for disc_prob in range(MAX_DISCRETE_PROB - 1, MIN_DISCRETE_PROB - 1, -1):
        p_values[disc_prob - MIN_DISCRETE_PROB] = (p_values[disc_prob + 1 - MIN_DISCRETE_PROB]
                                                   + to_prob(disc_prob) * path_counts[disc_prob - MIN_DISCRETE_PROB])

    # Subtract the contribution of one path with each discrete probability
    contributions = np.array([to_prob(DiscreteProb(k))
                              for k in range(MIN_DISCRETE_PROB, MAX_DISCRETE_PROB + 1)])
    p_values -= contributions

    # Prevent possible miscalculation caused by discretization
    p_values = np.clip(p_values, 0.0, 1.0)

    return p_values


def compare_p_values(hmm: HMM) -> Tuple[float, float, Dict[float, int]]:
    estimator = PValueEstimator(hmm)
    optimal_p_values = estimator._p_values

    true_p_values = compute_naive_p_values(hmm)

    # Calculate absolute differences
    absolute_differences = np.abs(optimal_p_values - true_p_values)

    max_absolute_difference = np.max(absolute_differences)
    mean_absolute_difference = np.mean(absolute_differences)

    # Create a dictionary of differences and their counts
    diff_dict = {}
    for diff in absolute_differences:
        rounded_diff = round(diff, 6)  # Round to 6 decimal places
        if rounded_diff in diff_dict:
            diff_dict[rounded_diff] += 1
        else:
            diff_dict[rounded_diff] = 1

    return max_absolute_difference, mean_absolute_difference, diff_dict


def create_test_hmm() -> HMM:
    transitions = [
        [(StateIdx(1), LogProb(-1.1)), (StateIdx(2), LogProb(-0.7)), (StateIdx(3), LogProb(-1.6))],
        [(StateIdx(4), LogProb(-0.4)), (StateIdx(5), LogProb(-1.0))],
        [(StateIdx(4), LogProb(-1.2)), (StateIdx(5), LogProb(-0.8)), (StateIdx(6), LogProb(-1.3))],
        [(StateIdx(5), LogProb(-0.3)), (StateIdx(6), LogProb(-0.9))],
        [(StateIdx(6), LogProb(-0.5)), (StateIdx(7), LogProb(-0.9))],
        [(StateIdx(7), LogProb(-0.2))],
        [(StateIdx(7), LogProb(-0.4))],
        []
    ]

    emissions = [
        [],
        [LogProb(-0.3), LogProb(-1.2)],
        [LogProb(-0.7), LogProb(-1.1), LogProb(-0.9)],
        [LogProb(-0.2), LogProb(-1.6)],
        [LogProb(-0.5)],
        [LogProb(-0.4), LogProb(-0.8)],
        [LogProb(-0.6), LogProb(-1.0), LogProb(-1.3)],
        []
    ]

    return HMM(transitions=transitions, emissions=emissions)


def test_p_value_estimator():
    test_hmm = create_test_hmm()

    max_abs_diff, mean_abs_diff, diff_dict = compare_p_values(test_hmm)

    print("\nDistribution of absolute differences:")
    print("Difference | Count")
    print("-" * 20)

    # Sort the differences for clearer output
    for diff, count in sorted(diff_dict.items()):
        print(f"{diff:10.6f} | {count}")
    print("-" * 20)
    print(f"  Total    | {MAX_DISCRETE_PROB - MIN_DISCRETE_PROB + 1}")

    print("-" * 20)
    print(f"Maximum absolute difference: {max_abs_diff:.6f}")
    print(f"Mean absolute difference: {mean_abs_diff:.6f}")


if __name__ == "__main__":
    test_p_value_estimator()