#include <cmath>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <iostream>
#include "data_types.h"
#include "p_values/p_values.h"

using namespace std;

using DiscreteLogOdds = int;
using LogProb = double; // log probability -- should be nonpositive
using LogOdds = double;
using Prob = double;  // probability -- should be between 0 and 1

const DiscreteLogOdds MIN_DISCRETE_LOG_ODDS = 0;
const DiscreteLogOdds MAX_DISCRETE_LOG_ODDS = 50000;
const LogProb LOG_ODDS_MIN_VALUE = -20.0;
const LogProb LOG_ODDS_MAX_VALUE = 50.0;

DiscreteLogOdds to_discrete_log_odds(LogOdds log_odds) {
    if (log_odds < LOG_ODDS_MIN_VALUE) {
        return MIN_DISCRETE_LOG_ODDS;
    }
    if (log_odds > LOG_ODDS_MAX_VALUE) {
        return MAX_DISCRETE_LOG_ODDS;
    }

    return static_cast<DiscreteLogOdds>(round(
            MIN_DISCRETE_LOG_ODDS + (log_odds - LOG_ODDS_MIN_VALUE) /
                                    (LOG_ODDS_MAX_VALUE - LOG_ODDS_MIN_VALUE) * (MAX_DISCRETE_LOG_ODDS - MIN_DISCRETE_LOG_ODDS)
    ));
}

LogOdds to_log_odds(DiscreteLogOdds disc_log_odds) {
    if (disc_log_odds == MIN_DISCRETE_LOG_ODDS) {
        return -INFINITY;
    } else {
        return LOG_ODDS_MIN_VALUE +
               (static_cast<double>(disc_log_odds - MIN_DISCRETE_LOG_ODDS) / (MAX_DISCRETE_LOG_ODDS - MIN_DISCRETE_LOG_ODDS)) *
               (LOG_ODDS_MAX_VALUE - LOG_ODDS_MIN_VALUE);
    }
}

vector<double> compute_hmm_p_values(const HMM& hmm_lo, const HMM& hmm_lp) {

    int num_states = hmm_lp.transitions.size();
    int initial_state = 0;
    int final_state = num_states - 1;

    vector<vector<Prob>> dp(num_states, vector<Prob>(MAX_DISCRETE_LOG_ODDS + 1, 0));
    dp[initial_state][to_discrete_log_odds(0.0)] = 1.0;

    for (int state = 0; state < num_states; ++state) {
        for (int disc_log_odds = MAX_DISCRETE_LOG_ODDS; disc_log_odds > MIN_DISCRETE_LOG_ODDS; --disc_log_odds) {
            if (dp[state][disc_log_odds] == 0.0) continue;

            Prob current_prob = dp[state][disc_log_odds];
            LogOdds current_lo = to_log_odds(disc_log_odds);

            for (int edge_idx = 0; edge_idx < hmm_lp.transitions[state].size(); edge_idx++){
                const auto& [next_state, transition_lp] = hmm_lp.transitions[state][edge_idx];
                const auto transition_lo = hmm_lo.transitions[state][edge_idx].second;
                if (!hmm_lp.emissions[next_state].empty()) {
                    for(int emission_idx = 0; emission_idx < hmm_lp.emissions[next_state].size(); emission_idx++){
                        LogProb emission_lp = hmm_lp.emissions[next_state][emission_idx];
                        LogOdds emission_lo = hmm_lo.emissions[next_state][emission_idx];

                        LogOdds new_lo = current_lo + transition_lo + emission_lo;
                        DiscreteLogOdds disc_new_lo = to_discrete_log_odds(new_lo);
                        Prob new_prob = current_prob * exp(emission_lp) * exp(transition_lp);
                        dp[next_state][disc_new_lo] += new_prob;
                    }
                } else {
                    LogProb new_lo = current_lo + transition_lo;
                    DiscreteLogOdds disc_new_lo = to_discrete_log_odds(new_lo);
                    Prob number_of_paths = current_prob * exp(transition_lp);
                    dp[next_state][disc_new_lo] += number_of_paths;
                }
            }
        }
    }

    vector<Prob> p_values(MAX_DISCRETE_LOG_ODDS - MIN_DISCRETE_LOG_ODDS + 1, 0.0);

    p_values[MAX_DISCRETE_LOG_ODDS] = dp[final_state][MAX_DISCRETE_LOG_ODDS];

    for (int disc_log_odds = MAX_DISCRETE_LOG_ODDS - 1; disc_log_odds >= MIN_DISCRETE_LOG_ODDS; --disc_log_odds) {
        p_values[disc_log_odds] = p_values[disc_log_odds + 1] + dp[final_state][disc_log_odds];
    }

    for (double& val : p_values) {
        val = min(max(val, 0.0), 1.0);
    }

    return p_values;
}
