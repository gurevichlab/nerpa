#pragma once

#include <vector>
#include <string>
#include <tuple>

// ========== Basic Type Aliases ==========

using MonCode    = int;
using rBAN_idx   = int;
using LogProb    = float;
using StateIdx   = int;
using NRP_ID     = std::string;
using BGC_Info   = std::tuple<std::string, int, int, int>;  // (genome_id, contig_idx, bgc_idx, variant_idx)

// ========== Data Structures ==========

// For each linearization, store sequences of monomers (first) and rBAN indices (second).
using NRP_Linearization = std::pair<std::vector<MonCode>, std::vector<rBAN_idx>>;

struct NRP_Linearizations {
    std::vector<NRP_Linearization> non_iterative;
    std::vector<std::vector<std::vector<NRP_Linearization>>> iterative;
};

// HMM data structure
struct HMM {
    // transitions[i]: vector of (next_state, log_prob)
    std::vector<std::vector<std::pair<StateIdx, LogProb>>> transitions;
    // emissions[i][j]: log probability of emitting monomer j from state i
    std::vector<std::vector<LogProb>> emissions;
    // nearest_module_start_state[i]: the nearest module start for state i
    std::vector<StateIdx> nearest_module_start_state;
};

struct MatchingConfig {
    int max_num_matches_per_bgc;
    int max_num_matches_per_nrp;
    int max_num_matches;
};

// ========== MatchInfo Structures ==========

struct MatchInfo {
    LogProb score;
    BGC_Info bgc_info;
    NRP_ID nrp_id;
    std::vector<NRP_Linearization> linearizations;
    std::vector<std::vector<StateIdx>> optimal_paths;  // one path per linearization
};

struct MatchInfoLight {
    LogProb score;
    BGC_Info bgc_info;
    NRP_ID nrp_id;
    std::vector<NRP_Linearization*> linearizations;
};

// ========== Forward Declarations (Optional) ==========
// If you need them, for example:
// LogProb get_hmm_score(...);
// std::pair<LogProb, std::vector<StateIdx>> get_opt_hmm_path(...);
