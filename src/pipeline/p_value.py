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


def discretization_error(hmm: HMM) -> Prob:
    # Error of discretization of one transition or emission
    single_error = 0.5 / (MAX_DISCRETE_PROB - MIN_DISCRETE_PROB)

    # Longest path contains less than number of vertices transitions
    transition_error = single_error * (len(hmm.transitions) - 1)

    # Longest path also contains emissions
    emission_error = 0.0
    for emission in hmm.emissions:
        if emission:
            emission_error += single_error

    return transition_error + emission_error


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


def compare_p_values(hmm: HMM) -> List[Tuple[DiscreteProb, LogProb, Prob, Prob, Prob, bool]]:
    # Get the optimal p-values
    estimator = PValueEstimator(hmm)
    optimal_p_values = estimator._p_values

    # Calculate the discretization error once
    error = discretization_error(hmm)

    # Get all paths and their log probabilities
    path_log_probs = get_all_paths_log_probs(hmm)
    path_probs = [np.exp(float(lp)) for lp in path_log_probs]

    # Results container
    results = []

    # Check each discrete probability
    for disc in range(MIN_DISCRETE_PROB, MAX_DISCRETE_PROB + 1):
        disc_prob = DiscreteProb(disc)

        # Get the optimal p-value (convert from numpy float to Prob)
        optimal_p_value = Prob(float(optimal_p_values[disc]))

        # Convert discrete probability to log probability
        threshold_log_prob = to_log_prob(disc_prob)

        # Calculate the lower and upper bounds in log probability space
        lower_bound_log_prob = LogProb(float(threshold_log_prob) - float(error))
        upper_bound_log_prob = LogProb(float(threshold_log_prob) + float(error))

        # Convert to actual probabilities
        lower_threshold_prob = np.exp(float(lower_bound_log_prob))
        upper_threshold_prob = np.exp(float(upper_bound_log_prob))

        # Calculate p-values for bounds by summing probabilities above threshold
        lower_bound_p_value = sum(p for p in path_probs if p >= lower_threshold_prob)
        upper_bound_p_value = sum(p for p in path_probs if p >= upper_threshold_prob)

        # Subtract one path contribution if there's an exact match
        exact_lower_matches = [p for p in path_probs if abs(p - lower_threshold_prob) < 1e-10]
        if exact_lower_matches:
            lower_bound_p_value -= exact_lower_matches[0]

        exact_upper_matches = [p for p in path_probs if abs(p - upper_threshold_prob) < 1e-10]
        if exact_upper_matches:
            upper_bound_p_value -= exact_upper_matches[0]

        # Ensure p-values are in valid range and convert to Prob type
        lower_bound_p_value = Prob(max(0.0, min(1.0, lower_bound_p_value)))
        upper_bound_p_value = Prob(max(0.0, min(1.0, upper_bound_p_value)))

        # Check if the optimal p-value is within bounds
        is_within_bounds = (float(lower_bound_p_value) <= float(optimal_p_value) <= float(upper_bound_p_value))

        # Store the results
        results.append((disc_prob, threshold_log_prob, optimal_p_value,
                        lower_bound_p_value, upper_bound_p_value, is_within_bounds))

    return results


def test_p_value_estimator():
    hmm = create_test_hmm()

    error = discretization_error(hmm)

    comparison_results = compare_p_values(hmm)

    # Print header
    print(f"Discretization error: {error}")
    print("\nComparing p-values for each threshold:")
    print("-" * 80)
    print(
        f"{'Disc Prob':<10} {'Log Prob':<15} {'Optimal p-value':<20} {'Lower Bound':<15} {'Upper Bound':<15} {'Within Bounds'}")
    print("-" * 80)

    # Print each result
    for disc_prob, log_prob, optimal_p, lower_bound, upper_bound, is_within in comparison_results:
        print(
            f"{disc_prob:<10} {float(log_prob):<15.6f} {float(optimal_p):<20.6e} {float(lower_bound):<15.6e} {float(upper_bound):<15.6e} {'Yes' if is_within else 'No'}")

    # Check if all values are within bounds
    all_within_bounds = all(result[5] for result in comparison_results)

    return all_within_bounds


if __name__ == "__main__":
    test_p_value_estimator()