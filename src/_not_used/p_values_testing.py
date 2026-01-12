

class PValueComparisonResult(NamedTuple):
    disc_prob: DiscreteLogOdds  # Discretized probability
    threshold_log_prob: LogProb  # Lof probability threshold
    true_p_value: Prob  # True p-value
    approx_p_value: Prob  # Approximate p-value
    lower_bound_p_value: Prob  # Lower bound of p-value
    upper_bound_p_value: Prob  # Upper bound of p-value
    is_within_bounds: bool  # Whether the approximate p-value is within the bounds

def discretization_error(hmm: HMM) -> LogProb:
    # Error of discretization of one transition or emission
    single_error = 0.5 * (LOG_ODDS_MAX_VALUE - LOG_ODDS_MIN_VALUE + 1) / (MAX_DISCRETE_LOG_ODDS - MIN_DISCRETE_LOG_ODDS + 1)

    # Longest path contains less than number of vertices transitions
    transition_error = single_error * (len(hmm.transitions) - 1)

    # Longest path also contains emissions
    emission_error = 0.0
    for state_emissions in hmm.emissions:
        if state_emissions:
            emission_error += single_error

    return transition_error + emission_error



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
        [(StateIdx(1), LogProb(-1.0)), (StateIdx(2), LogProb(-1.2)), (StateIdx(3), LogProb(-1.5))],
        [(StateIdx(4), LogProb(-0.8)), (StateIdx(5), LogProb(-1.1)), (StateIdx(6), LogProb(-1.3))],
        [(StateIdx(4), LogProb(-1.0)), (StateIdx(5), LogProb(-0.9)), (StateIdx(6), LogProb(-1.2))],
        [(StateIdx(4), LogProb(-1.3)), (StateIdx(5), LogProb(-1.0)), (StateIdx(6), LogProb(-0.7))],
        [(StateIdx(7), LogProb(-0.6)), (StateIdx(8), LogProb(-1.2)), (StateIdx(9), LogProb(-1.0))],
        [(StateIdx(7), LogProb(-0.9)), (StateIdx(8), LogProb(-0.8)), (StateIdx(9), LogProb(-1.1))],
        [(StateIdx(7), LogProb(-1.2)), (StateIdx(8), LogProb(-0.7)), (StateIdx(9), LogProb(-0.9))],
        [(StateIdx(10), LogProb(-0.3))],
        [(StateIdx(10), LogProb(-0.4))],
        [(StateIdx(10), LogProb(-0.5))],
        []
    ]

    emissions = [
        [],
        [LogProb(-0.5), LogProb(-1.2)],
        [LogProb(-0.6), LogProb(-1.0), LogProb(-1.4)],
        [LogProb(-0.3), LogProb(-1.1)],
        [LogProb(-0.7)],
        [LogProb(-0.4), LogProb(-0.9)],
        [LogProb(-0.6), LogProb(-1.3), LogProb(-0.8)],
        [LogProb(-0.2), LogProb(-1.5)],
        [LogProb(-0.5), LogProb(-1.0)],
        [LogProb(-0.8)],
        []
    ]

    return HMM(transitions=transitions, emissions=emissions)


def prepare_p_value_comparison(hmm: HMM) -> Tuple[PValueEstimator, LogProb, List[float]]:
    """
    Computes p-values, creates the p-value estimator, and calculates the discretization error.
    """
    # Create the p-value estimator
    p_value_estimator = PValueEstimator(hmm)

    # Calculate the discretization error
    log_error = discretization_error(hmm)

    # Get all paths and their log probabilities
    path_log_probs = get_all_paths_log_probs(hmm)
    path_probs = [np.exp(float(lp)) for lp in path_log_probs]

    return p_value_estimator, log_error, path_probs


def check_p_value_bounds(p_value_estimator: PValueEstimator, log_error: LogProb, path_probs: List[float]) -> List[PValueComparisonResult]:
    """
    Checks if the approximate p-values are within bounds.
    """
    error = log_error
    results = []

    # Check each discrete probability
    for disc in range(MIN_DISCRETE_LOG_ODDS, MAX_DISCRETE_LOG_ODDS + 1):
        disc_prob = DiscreteLogOdds(disc)

        # Convert discrete probability to log probability
        threshold_log_prob = to_log_odds(disc_prob)

        # Get the optimal p-value using the estimator
        optimal_p_value = p_value_estimator(threshold_log_prob)

        # Calculate the lower and upper bounds in log probability space
        threshold_lp_lb = LogProb(float(threshold_log_prob - log_error))
        threshold_lp_ub = LogProb(float(threshold_log_prob + log_error))

        lower_threshold_prob = np.exp(float(threshold_lp_lb))
        upper_threshold_prob = np.exp(float(threshold_lp_ub))

        # Reverse the bounds since p-value(threshold) is decreasing function
        upper_bound_p_value = sum(p for p in path_probs if p >= lower_threshold_prob)
        lower_bound_p_value = sum(p for p in path_probs if p >= upper_threshold_prob)

        log_prob = to_log_odds(disc_prob)
        actual_prob = np.exp(float(log_prob))

        if error > actual_prob - 1e-5:
            upper_error = float('inf')
        else:
            upper_error = error / (actual_prob - error)
        lower_error = error / (actual_prob + error)

        # Account the second type of error caused by crossing the threshold by probabilities of paths
        upper_bound_p_value += upper_error
        lower_bound_p_value -= lower_error

        # "<=" is used to include the case when p-value is zero or one
        is_within_bounds = (float(lower_bound_p_value) <= float(optimal_p_value) <= float(upper_bound_p_value))

        true_p_value = sum(p for p in path_probs if p >= actual_prob)

        # Store the results
        results.append(PValueComparisonResult(
            disc_prob, threshold_log_prob, true_p_value, optimal_p_value,
            lower_bound_p_value, upper_bound_p_value, is_within_bounds
        ))

    return results


def test_p_value_estimator():
    hmm = create_test_hmm()

    # Prepare p-value comparison components
    p_value_estimator, log_error, path_probs = prepare_p_value_comparison(hmm)

    # Check if p-values are within bounds
    comparison_results = check_p_value_bounds(p_value_estimator, log_error, path_probs)

    # Print header
    with open('p_values_test.txt', 'w') as f:
        print((f"Discretization error: {log_error}\n"
               + "-" * 100
               + f"\n{'Disc Prob':<10} {'Log Prob':<15} {'True p-value':<15} {'Approx p-value':<20}"
               + f" {'Lower Bound':<15} {'Upper Bound':<15} {'Within Bounds'}\n"
               + "-" * 100),
              file=f)

        # Print each result
        step = int(len(comparison_results) / 100)
        for i in range(100):
            result = comparison_results[step * i]
            print((
                f"{result.disc_prob:<10} {float(result.threshold_log_prob):<15.6f} "
                f"{float(result.true_p_value):<15.6f} "
                f"{float(result.approx_p_value):<20.6f} {float(result.lower_bound_p_value):<15.6f} "
                f"{float(result.upper_bound_p_value):<15.6f} {'Yes' if result.is_within_bounds else 'No'}"
            ), file=f)

    # Check if all values are within bounds
    all_within_bounds = all(result.is_within_bounds for result in comparison_results)

    return all_within_bounds


if __name__ == "__main__":
    test_p_value_estimator()
