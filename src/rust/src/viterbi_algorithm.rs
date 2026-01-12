use crate::data_types::*;

pub fn get_hmm_score(hmm: &HMM, nrp_monomers: &[MonCode]) -> LogProb {
    let seq_len: usize   = nrp_monomers.len();                 // number of observed symbols
    let num_states: usize = hmm.num_states();
    let initial_state: usize = 0;
    let final_state: usize = num_states - 1;

    // dp[state][k] with k = #symbols consumed (0..=seq_len)
    let mut dp: Vec<Vec<LogProb>> = vec![vec![NEG_INFINITY; seq_len + 1]; num_states];
    dp[initial_state][0] = ZERO;

    //let alpha = hmm.alphabet();
    //let obs_usize: Vec<usize> = nrp_monomers.iter().map(|&m| m as usize).collect();

    for num_symb_processed in 0..=seq_len {
        /* If we ignore gene/module iteration edges, states are sorted topologically.
       After an iteration edge, however, the state index is smaller than the previous one.
       Between any two iterations, a symbol should be emitted (otherwise iteration is useless).
       So I iterate over all states, and if an iteration edge is applied,
       I iterate one more time over the same range. */

        let mut iteration_applied = false;
        for _dp_run in 0..2 {
            if _dp_run == 1 && !iteration_applied { break; }  // no need for second run if no iteration edge was applied

            for state_idx in 0..num_states {
                if (num_symb_processed == seq_len
                    && !hmm.emissions[state_idx].is_empty()) { continue; }

                let log_prob = dp[state_idx][num_symb_processed];  // best score for this state and #symbols processed
                if log_prob == NEG_INFINITY { continue; }  // the dp state is unreachable

                // Process all outgoing transitions from the state
                for &(next_state_idx, edge_log_prob) in &hmm.transitions[state_idx] {
                    let mut new_num_symbols_processed = num_symb_processed;
                    let mut new_log_prob = log_prob + edge_log_prob;
                    if !hmm.emissions[state_idx].is_empty() {
                        let emitted_symbol = nrp_monomers[num_symb_processed];
                        new_log_prob += hmm.emissions[state_idx][emitted_symbol];
                        new_num_symbols_processed += 1;
                    }
                    if new_log_prob > dp[next_state_idx][new_num_symbols_processed] {
                        dp[next_state_idx][new_num_symbols_processed] = new_log_prob;
                        if next_state_idx < state_idx && new_num_symbols_processed == num_symb_processed {
                            iteration_applied = true; // match back-edge trigger
                        }
                    }
                }
            }
        }
    }

    dp[final_state][seq_len]
}

pub fn get_hmm_score_with_path(hmm: &HMM, nrp_monomers: &[MonCode]) -> (LogProb, Vec<StateIdx>) {
    let seq_len = nrp_monomers.len();
    let num_states = hmm.num_states();
    let initial_state = 0;
    let final_state = num_states - 1;

    let mut dp: Vec<Vec<LogProb>> = vec![vec![NEG_INFINITY; seq_len + 1]; num_states];
    let mut prev: Vec<Vec<Option<StateIdx>>>   = vec![vec![None; seq_len + 1]; num_states];
    dp[initial_state][0] = ZERO;

    for num_symb_processed in 0..=seq_len {
        let mut iteration_applied = false;
        for _dp_run in 0..2 {
            if _dp_run == 1 && !iteration_applied { break; }
            for state_idx in 0..num_states {
                if num_symb_processed == seq_len && !hmm.emissions[state_idx].is_empty() { continue; }
                let log_prob = dp[state_idx][num_symb_processed];
                if log_prob == NEG_INFINITY { continue; }
                for &(next_state_idx, edge_log_prob) in &hmm.transitions[state_idx] {
                    let mut new_num_symbols_processed = num_symb_processed;
                    let mut new_log_prob = log_prob + edge_log_prob;
                    if !hmm.emissions[state_idx].is_empty() {
                        let emitted_symbol = nrp_monomers[num_symb_processed];
                        new_log_prob += hmm.emissions[state_idx][emitted_symbol];
                        new_num_symbols_processed += 1;
                    }
                    if new_log_prob > dp[next_state_idx][new_num_symbols_processed] {
                        dp[next_state_idx][new_num_symbols_processed] = new_log_prob;
                        prev[next_state_idx][new_num_symbols_processed] = Some(state_idx);
                        if next_state_idx < state_idx && new_num_symbols_processed == num_symb_processed {
                            iteration_applied = true;
                        }
                    }
                }
            }
        }
    }

    let score = dp[final_state][seq_len];

    // Reconstruct path: (final_state, seq_len) -> (initial_state, 0)
    let mut path: Vec<StateIdx> = Vec::new();
    let mut cur_state: StateIdx = final_state;
    let mut cur_num_processed_symbols = seq_len;
    path.push(cur_state);
    while !(cur_state == initial_state && cur_num_processed_symbols == 0) {
        let prev_state = prev[cur_state][cur_num_processed_symbols]
            .expect("previous state should always exist while backtracking");
        if !hmm.emissions[prev_state].is_empty() { cur_num_processed_symbols -= 1; } // If the previous state emits, decrement the symbol count.
        cur_state = prev_state;
        path.push(cur_state);
    }
    path.reverse();

    // Validate the endpoints (like C++ asserts)
    debug_assert!(!path.is_empty());
    debug_assert_eq!(*path.first().unwrap(), initial_state, "the path endpoints are wrong!!!");
    debug_assert_eq!(*path.last().unwrap(),  final_state,   "the path endpoints are wrong!!!");

    (score, path)
}