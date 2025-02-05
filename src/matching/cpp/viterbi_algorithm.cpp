#include "data_types.h"
using namespace std;


// checkpoint (state_idx, num_symbols) means that by reaching state_idx, num_symbols symbols should have been emitted
LogProb get_hmm_score(const HMM& hmm,
                      const vector<MonCode>& nrp_monomers,
                      const vector<pair<StateIdx, int>>& checkpoints) {
    int seq_len = observed_sequence.size();         // length of the sequence
    int num_states = hmm.adj_list.size();                // number of states

    // dp[u][i] = maximum log probability ending in state u after processing i symbols
    vector<vector<double>> dp(m, vector<LogProb>(seq_len + 1, -numeric_limits<double>::infinity()));

    // Starting condition: at first checkpoint, the log probability is 0.
    dp[checkpoints[0].first][checkpoints[0].second] = 0.0;

    // Each element in the queue is a triple: (state, num symbols processed, current log probability)
    using Triple = tuple<StateIdx , int, LogProb>;
    queue<Triple> q;
    q.emplace(checkpoints[0].first, checkpoints[0].second, 0.0);

    for(const auto [checkpnt_state_idx, checkpnt_num_symb_processed] : checkpoints) {
        // on each iteration, we find the best path until the next checkpoint
        while(!q.empty()) {
            const auto [edge_from, num_symb_processed, log_prob] = q.front();
            q.pop();

            // If we've processed all symbols and the current state should emit a symbol,
            if (num_symb_processed == checkpnt_num_symb_processed && !hmm.emission_log_probs[edge_from].empty()) {
                continue;
            }

            // Iterate over all outgoing edges from state edge_from.
            for (const auto [edge_to, edge_log_prob] : hmm.adj_list[edge_from]) {
                if hmm.nearest_module_start_state[edge_to] > hmm.nearest_module_start_state[checkpnt_state_idx] {
                    continue;  // next checkpoint is unreachable from edge_to
                }
                auto [new_num_symb_processed, new_log_prob] = num_symb_processed, log_prob + edge_log_prob;

                // If the state emits a symbol, update the log probability with the emission probability.
                if (!hmm.emission_log_probs[edge_from].empty()) {
                    // It is assumed that num_symb_processed < n here because otherwise we would have continued above.
                    new_log_prob += hmm.emission_log_probs[edge_from][observed_sequence[num_symb_processed]];
                    new_num_symb_processed++;
                }

                // Only update if this new path gives a higher log probability.
                if (new_log_prob > dp[edge_to][new_num_symb_processednum_symbols]) {
                    dp[edge_to][new_num_symb_processed] = new_log_prob;
                    q.emplace(edge_to, new_num_symb_processed, new_log_prob);
                }
            }
        }
        // start new iteration from the next checkpoint
        q.emplace(checkpnt_state_idx,
                  checkpnt_num_symb_processed,
                  dp[checkpnt_state_idx][checkpnt_num_symb_processed]);
    }

    }

    return dp[checkpoints.back().first][checkpoints.back().second];
}