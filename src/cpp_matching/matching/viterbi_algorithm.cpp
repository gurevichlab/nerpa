#include <cassert>
#include "data_types.h"
#include <vector>
#include <limits>
#include <utility>
using namespace std;

// Computes the maximum probability path through an HMM given an input sequence.
LogProb get_hmm_score(const HMM& hmm,
                      const vector<MonCode>& nrp_monomers,
                      const vector<pair<StateIdx, int>>& checkpoints) {
    int seq_len = nrp_monomers.size();       // Number of observed symbols
    int num_states = hmm.transitions.size(); // Number of states in the HMM

    // dp[state][num_symbols] = max log probability reaching `state` after processing `num_symbols` monomers
    vector<vector<LogProb>> dp(num_states, vector<LogProb>(seq_len + 1,
                                                           -numeric_limits<LogProb>::infinity()));

    // Initialize at first checkpoint
    dp[checkpoints[0].first][checkpoints[0].second] = 0.0;

    // Iterate through all checkpoints except the last (final state is already covered)
    for (size_t checkpoint_idx = 0; checkpoint_idx < checkpoints.size() - 1; ++checkpoint_idx) {
        auto [checkpnt_state_idx, checkpnt_num_symb_processed] = checkpoints[checkpoint_idx];
        auto [next_checkpnt_state_idx, next_checkpnt_num_symb_processed] = checkpoints[checkpoint_idx + 1];

        for (int num_symb_processed = checkpnt_num_symb_processed;
             num_symb_processed <= next_checkpnt_num_symb_processed; ++num_symb_processed) {

            /* If we ignore gene/module iteration edges, states are sorted topologically.
            After an iteration edge, however, the state index is less than the previous one.
            Between any two iterations, a symbol should be emitted (otherwise iteration is useless).
            So I iterate over all states in the range [cur_state, next_state) and if an iteration edge is applied,
            I iterate one more time over the same range. This is done at most once. */

            bool iteration_applied = false;
            for(int _dp_run = 0; _dp_run < 2; _dp_run++) {
                if (_dp_run == 1 and not iteration_applied) break;

                for (StateIdx state_idx = checkpnt_state_idx; state_idx < next_checkpnt_state_idx; ++state_idx) {
                    if (num_symb_processed == next_checkpnt_num_symb_processed and
                        !hmm.emissions[state_idx].empty()) {
                        continue;
                    }

                    LogProb log_prob = dp[state_idx][num_symb_processed];
                    if (log_prob == -numeric_limits<LogProb>::infinity()) continue; // Skip uninitialized states

                    for (const auto &[edge_to, edge_log_prob]: hmm.transitions[state_idx]) {
                        int new_num_symb_processed = num_symb_processed;
                        LogProb new_log_prob = log_prob + edge_log_prob;

                        if (!hmm.emissions[state_idx].empty()) { // If the state emits a symbol
                            new_log_prob += hmm.emissions[state_idx][nrp_monomers[num_symb_processed]];
                            new_num_symb_processed++;
                        }

                        if (new_log_prob > dp[edge_to][new_num_symb_processed]) {
                            dp[edge_to][new_num_symb_processed] = new_log_prob;
                            if (edge_to < state_idx and new_num_symb_processed == num_symb_processed) {
                                iteration_applied = true;
                            }
                        }
                    }
                }
            }
        }
    }

    // Return the probability of reaching the final state after processing all symbols
    return dp[checkpoints.back().first][checkpoints.back().second];
}

pair<LogProb, vector<StateIdx>>
get_opt_hmm_path(const HMM& hmm,
                 StateIdx initial_state,
                 StateIdx final_state,
                 const vector<MonCode>& nrp_monomers) {
    int seq_len = static_cast<int>(nrp_monomers.size());  // number of observed symbols
    int num_states = static_cast<int>(hmm.transitions.size());  // number of states


    // dp[u][i] = maximum log probability ending in state u after processing i symbols
    vector<vector<LogProb>> dp(num_states,
                               vector<LogProb>(seq_len + 1, -numeric_limits<LogProb>::infinity()));
    // prev[u][i] = previous state for state u when i symbols have been processed (for path reconstruction)
    vector<vector<StateIdx>> prev(num_states, vector<StateIdx>(seq_len + 1, -1));

    // Starting condition.
    dp[initial_state][0] = 0.0;

    // Iterate through all checkpoints except the last (final state is already covered)
    for (int num_symb_processed = 0; num_symb_processed <= seq_len; ++num_symb_processed) {
        bool iteration_applied = false;  // see get_hmm_score for explanation
        for (int _dp_run = 0; _dp_run < 2; _dp_run++) {
            if (_dp_run == 1 and not iteration_applied) break;

            for (StateIdx state_idx = 0; state_idx < num_states; ++state_idx) {

                if (num_symb_processed == seq_len and !hmm.emissions[state_idx].empty()) continue;

                LogProb log_prob = dp[state_idx][num_symb_processed];
                if (log_prob == -numeric_limits<LogProb>::infinity()) continue; // Skip uninitialized states

                // Process all outgoing transitions from the state
                for (const auto &[next_state_idx, edge_log_prob]: hmm.transitions[state_idx]) {
                    int new_num_symbols = num_symb_processed;
                    LogProb new_log_prob = log_prob + edge_log_prob;
                    // If the current state emits a symbol, incorporate its emission probability.
                    if (!hmm.emissions[state_idx].empty()) {
                        new_log_prob += hmm.emissions[state_idx][nrp_monomers[num_symb_processed]];
                        new_num_symbols++;
                    }
                    // Update if this new path gives a higher log probability.
                    if (new_log_prob > dp[next_state_idx][new_num_symbols]) {
                        dp[next_state_idx][new_num_symbols] = new_log_prob;
                        prev[next_state_idx][new_num_symbols] = state_idx;
                        if (next_state_idx < state_idx and new_num_symbols == num_symb_processed) {
                            iteration_applied = true;
                        }
                    }
                }
            }
        }
    }

    // --- Reconstruct the optimal path ---
    StateIdx u = final_state;
    int k = seq_len;  // number of symbols processed when reaching final_state
    vector<StateIdx> path{final_state};
    while (u != initial_state) {
        StateIdx prev_state = prev[u][k];
        // Ensure that a previous state exists.
        assert(prev_state != -1 && "prev_state not defined!!!");
        // If the previous state emits a symbol, decrement the symbol count.
        if (!hmm.emissions[prev_state].empty()) k--;
        u = prev_state;
        path.push_back(u);
    }
    // Reverse the path since it was constructed backwards.
    reverse(path.begin(), path.end());

    // Validate the endpoints.
    assert(!path.empty());
    assert(path.front() == initial_state && "the path endpoints are wrong!!!");
    assert(path.back() == final_state && "the path endpoints are wrong!!!");

    return { dp[final_state][seq_len], path };
}