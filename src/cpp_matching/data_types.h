#pragma once

#include <vector>
#include <string>
#include <tuple>
#include <functional>
#include <utility>

// ========== Basic Type Aliases ==========

using MonCode    = int;
using rBAN_idx   = int;
using LogProb    = double;
using StateIdx   = int;
using NRP_ID     = std::string;

// ========== Data Structures ==========

using BGC_ID = std::tuple<std::string, int, int>;  // (genome_id, contig_idx, bgc_idx)
// q: write accessors for BGC_ID
constexpr const auto& genome_id(const BGC_ID & t) { return std::get<0>(t); }
constexpr const auto& contig_idx(const BGC_ID & t) { return std::get<1>(t); }
constexpr const auto& bgc_idx(const BGC_ID & t) { return std::get<2>(t); }

using BGC_Variant_ID  = std::tuple<BGC_ID, int>;  // (bgc_id, variant_idx)
// q: write accessors for BGC_Variant_ID
constexpr const auto& bgc_id(const BGC_Variant_ID & t) { return std::get<0>(t); }
constexpr const auto& variant_idx(const BGC_Variant_ID & t) { return std::get<1>(t); }

// For each linearization, store sequences of monomers (first) and rBAN indices (second).
using NRP_Linearization = std::pair<std::vector<MonCode>, std::vector<rBAN_idx>>;
constexpr auto& mon_codes(NRP_Linearization & t) { return std::get<0>(t); }
constexpr auto& rban_idxs(NRP_Linearization & t) { return std::get<1>(t); }


struct NRP_Linearizations {
    std::vector<NRP_Linearization> non_iterative;
    std::vector<std::vector<std::vector<NRP_Linearization>>> iterative;
    double score_vs_avg_bgc;  // used for ranking matches
};

// HMM data structure
struct HMM {
    BGC_Variant_ID  bgc_variant_id;
    // transitions[i]: vector of (next_state, log_prob)
    std::vector<std::vector<std::pair<StateIdx, LogProb>>> transitions;
    // emissions[i][j]: log probability of emitting monomer j from state i
    std::vector<std::vector<LogProb>> emissions;
    std::vector<StateIdx> module_start_states;
    std::vector<StateIdx> module_match_states;

};

struct MatchingConfig {
    int max_num_matches_per_bgc;
    int max_num_matches_per_nrp;
    int max_num_matches;
};

// ========== MatchInfo Structures ==========

struct MatchInfo {
    LogProb score;
    BGC_Variant_ID bgc_variant_id;
    NRP_ID nrp_id;
    std::vector<NRP_Linearization> linearizations;
    std::vector<std::vector<StateIdx>> optimal_paths;  // one path per linearization
};

struct MatchInfoLight {
    LogProb raw_score;
    LogProb score;  // used for ranking matches
    BGC_Variant_ID bgc_variant_id;
    NRP_ID nrp_id;
    std::vector<const NRP_Linearization*> linearizations;
};


// ========== Hashing for std::tuple ==========
namespace std {
    template <typename... T>
    struct hash<std::tuple<T...>> {
        std::size_t operator()(const std::tuple<T...>& tup) const {
            return hash_tuple(tup, std::index_sequence_for<T...>{});
        }

    private:
        template <typename Tuple, std::size_t... Is>
        static std::size_t hash_tuple(const Tuple& tup, std::index_sequence<Is...>) {
            std::size_t seed = 0;
            (..., (seed ^= std::hash<std::tuple_element_t<Is, Tuple>>{}(std::get<Is>(tup)) + 0x9e3779b9 + (seed << 6) + (seed >> 2)));
            return seed;
        }
    };
}