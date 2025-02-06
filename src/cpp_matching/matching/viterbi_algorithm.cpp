#include "../data_types.h"
#include <queue>
#include <tuple>
#include <limits>
#include <cassert>
#include <algorithm>
using namespace std;


// checkpoint (state_idx, num_symbols) means that by reaching state_idx, num_symbols symbols should have been emitted
LogProb get_hmm_score(const HMM& hmm,
                      const vector<MonCode>& nrp_monomers,
                      const vector<pair<StateIdx, int>>& checkpoints) {
    int seq_len = static_cast<int>(nrp_monomers.size());  // length of the sequence
    int num_states = static_cast<int>(hmm.transitions.size());  // number of states

    // dp[u][i] = maximum log probability ending in state u after processing i symbols
    vector<vector<LogProb>> dp(num_states,
                               vector<LogProb>(seq_len + 1, -numeric_limits<LogProb>::infinity()));

    // Starting condition: at the first checkpoint, the log probability is 0.
    dp[checkpoints[0].first][checkpoints[0].second] = 0.0;

    queue<tuple<StateIdx, int, LogProb>> q;
    q.emplace(checkpoints[0].first, checkpoints[0].second, 0.0);

    // For each checkpoint, find the best path until that checkpoint.
    for (const auto& [checkpnt_state_idx, checkpnt_num_symb_processed] : checkpoints) {
        while (!q.empty()) {
            auto [edge_from, num_symb_processed, log_prob] = q.front();
            q.pop();

            // If we've processed exactly the checkpoint's number of symbols and the current state emits,
            // then skip this path.
            if (num_symb_processed == checkpnt_num_symb_processed &&
                !hmm.emissions[edge_from].empty()) {
                continue;
            }

            // Iterate over all outgoing transitions from state 'edge_from'.
            for (const auto& [edge_to, edge_log_prob] : hmm.transitions[edge_from]) {
                // If the next state's nearest module start is greater than that of the checkpoint state,
                // then the next checkpoint is unreachable from edge_to.
                if (hmm.nearest_module_start_state[edge_to] >
                    hmm.nearest_module_start_state[checkpnt_state_idx]) {
                    continue;
                }
                // Initialize new variables.
                int new_num_symb_processed = num_symb_processed;
                LogProb new_log_prob = log_prob + edge_log_prob;

                // If the current state emits a symbol, update with its emission probability.
                if (!hmm.emissions[edge_from].empty()) {
                    // Assumes that num_symb_processed < seq_len.
                    new_log_prob += hmm.emissions[edge_from][nrp_monomers[num_symb_processed]];
                    new_num_symb_processed++;
                }

                // Update if this new path gives a higher log probability.
                if (new_log_prob > dp[edge_to][new_num_symb_processed]) {
                    dp[edge_to][new_num_symb_processed] = new_log_prob;
                    q.emplace(edge_to, new_num_symb_processed, new_log_prob);
                }
            }
        }
        // Start a new iteration from the current checkpoint.
        q.emplace(checkpnt_state_idx, checkpnt_num_symb_processed,
                  dp[checkpnt_state_idx][checkpnt_num_symb_processed]);
    }

    return dp[checkpoints.back().first][checkpoints.back().second];
}

pair<LogProb, vector<StateIdx>>
get_opt_hmm_path(const HMM& hmm,
                 StateIdx initial_state,
                 StateIdx final_state,
                 const vector<MonCode>& nrp_monomers) {
    int seq_len = static_cast<int>(nrp_monomers.size());  // number of observed symbols
    int num_states = static_cast<int>(hmm.transitions.size());  // number of states

    // Each element in the queue is a tuple: (current state, num symbols processed, current log probability)
    queue<tuple<StateIdx, int, LogProb>> q;
    q.emplace(initial_state, 0, 0.0);

    // dp[u][i] = maximum log probability ending in state u after processing i symbols
    vector<vector<LogProb>> dp(num_states,
                               vector<LogProb>(seq_len + 1, -numeric_limits<LogProb>::infinity()));
    // prev[u][i] = previous state for state u when i symbols have been processed (for path reconstruction)
    vector<vector<StateIdx>> prev(num_states, vector<StateIdx>(seq_len + 1, -1));

    // Starting condition.
    dp[initial_state][0] = 0.0;

    while (!q.empty()) {
        auto [u, num_symb_processed, log_prob] = q.front();
        q.pop();

        // If we've processed all symbols and the current state emits, skip this path.
        if (num_symb_processed == seq_len && !hmm.emissions[u].empty()) {
            continue;
        }

        // Process all outgoing transitions from state 'u'.
        for (const auto& [v, edge_log_prob] : hmm.transitions[u]) {
            int new_num_symbols = num_symb_processed;
            LogProb new_log_prob = log_prob + edge_log_prob;
            // If the current state emits a symbol, incorporate its emission probability.
            if (!hmm.emissions[u].empty()) {
                new_log_prob += hmm.emissions[u][nrp_monomers[num_symb_processed]];
                new_num_symbols++;
            }
            // Update if this new path gives a higher log probability.
            if (new_num_symbols <= seq_len && new_log_prob > dp[v][new_num_symbols]) {
                dp[v][new_num_symbols] = new_log_prob;
                prev[v][new_num_symbols] = u;
                q.emplace(v, new_num_symbols, new_log_prob);
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