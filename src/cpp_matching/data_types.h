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

using BGC_Info   = std::tuple<std::string, int, int, int>;  // (genome_id, contig_idx, bgc_idx, variant_idx)
// Define accessors
constexpr auto& genome_id(BGC_Info& t) { return std::get<0>(t); }
constexpr auto& genome_id(const BGC_Info& t) { return std::get<0>(t); }
constexpr auto& contig_idx(BGC_Info& t) { return std::get<1>(t); }
constexpr auto& bgc_idx(BGC_Info& t) { return std::get<2>(t); }
constexpr auto& variant_idx(BGC_Info& t) { return std::get<3>(t); }

// For each linearization, store sequences of monomers (first) and rBAN indices (second).
using NRP_Linearization = std::pair<std::vector<MonCode>, std::vector<rBAN_idx>>;
constexpr auto& mon_codes(NRP_Linearization & t) { return std::get<0>(t); }
constexpr auto& rban_idxs(NRP_Linearization & t) { return std::get<1>(t); }


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
    BGC_Info bgc_info;
    NRP_ID nrp_id;
    std::vector<NRP_Linearization> linearizations;
    std::vector<std::vector<StateIdx>> optimal_paths;  // one path per linearization
};

struct MatchInfoLight {
    LogProb score;
    BGC_Info bgc_info;
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