#include "../data_types.h"
#include "matcher.h"
#include "viterbi_algorithm.h"
#include "../openmp_wrapper.h"
#include <vector>
#include <limits>
#include <utility>
using namespace std;


LogProb linearization_score(const HMM& hmm,
                            const NRP_Linearization& linearization){
    vector<pair<StateIdx, int>> checkpoints = { {0, 0},
                                                {(int) hmm.transitions.size() - 1, linearization.first.size()} };
    return get_hmm_score(hmm, linearization.first, checkpoints);
}


std::pair<LogProb , std::vector<const NRP_Linearization*>>
get_best_linearizations_for_nrp(
        const HMM& hmm,
        const NRP_Linearizations& nrp_linearizations)
{
    // Determine the best non-iterative linearization.
    LogProb best_non_iterative_score = -std::numeric_limits<double>::infinity();
    const NRP_Linearization* best_non_iterative_linearization;  // Default constructed; assumes non-empty nrp_linearizations.non_iterative

    for (const auto& non_iterative_linearization : nrp_linearizations.non_iterative) {
        LogProb score = linearization_score(hmm, non_iterative_linearization);
        if (score > best_non_iterative_score) {
            best_non_iterative_score = score;
            best_non_iterative_linearization = &non_iterative_linearization;
        }
    }

    // Determine the best iterative linearizations.
    double best_iterative_score = -std::numeric_limits<double>::infinity();
    std::vector<const NRP_Linearization*> best_iterative_linearizations;

    // Assuming nrp_linearizations.iterative is a vector of "groups",
    // where each group is a vector (of groups) of linearizations.
    for (const auto& linearizations_of_one_split : nrp_linearizations.iterative) {
        double split_score = 0;
        std::vector<const NRP_Linearization*> split_linearizations;
        for (const auto& group_linearizations : linearizations_of_one_split) {
            double group_score = -std::numeric_limits<double>::infinity();
            const NRP_Linearization* best_group_linearization;
            for (const auto& linearization : group_linearizations) {
                double score = linearization_score(hmm, linearization);
                if (score > group_score) {
                    group_score = score;
                    best_group_linearization = &linearization;
                }
            }
            split_linearizations.push_back(best_group_linearization);
            split_score += group_score;
        }
        if (split_score > best_iterative_score) {
            best_iterative_score = split_score;
            best_iterative_linearizations = split_linearizations;
        }
    }


    // Return the best overall result.
    if (best_non_iterative_score > best_iterative_score) {
        return { best_non_iterative_score, { best_non_iterative_linearization } };
    } else {
        return { best_iterative_score, best_iterative_linearizations };
    }
}

vector<MatchInfoLight>
get_best_matches_for_bgc_variant(const HMM& hmm,
                                 const BGC_Info& bgc_info,
                                 const vector<pair<NRP_Linearizations, NRP_ID>>& all_nrp_linearizations,
                                 const int max_num_matches_per_bgc)
{
    vector<MatchInfoLight> best_matches;
    for (const auto& [nrp_linearizations, nrp_id] : all_nrp_linearizations) {
        auto [score, linearizations] = get_best_linearizations_for_nrp(hmm, nrp_linearizations);
        best_matches.push_back({score, bgc_info, nrp_id, linearizations});
    }

    sort(best_matches.begin(), best_matches.end(),
         [](const auto& a, const auto& b) {
             return a.score > b.score;
         });

    // Truncate to at most max_num_matches_per_bgc entries.
    if (best_matches.size() > static_cast<size_t>(max_num_matches_per_bgc) and max_num_matches_per_bgc > 0) {
        best_matches.resize(max_num_matches_per_bgc);
    }
    return best_matches;
}


vector<MatchInfoLight> filter_and_sort_matches(const vector<MatchInfoLight>& matches,
                                               const MatchingConfig& config)
{
    {
        // Step 0: For each NRP, keep only the top variant per BGC
        unordered_map<NRP_ID, unordered_map<BGC_Info, MatchInfoLight>> unique_matches_map;
        for (const auto& match : matches) {
            // Normalize the BGC key by setting its variant_idx to 0.
            BGC_Info normalizedBGC = match.bgc_info;
            variant_idx(normalizedBGC) = 0;

            // Insert or update the best match for this (nrp_id, normalizedBGC) pair.
            auto& bgc_map = unique_matches_map[match.nrp_id];
            auto it = bgc_map.find(normalizedBGC);
            if (it == bgc_map.end() || match.score > it->second.score) {
                bgc_map[normalizedBGC] = match;
            }
        }

        // Collect the unique matches.
        vector<MatchInfoLight> unique_matches;
        for (const auto& nrp_pair : unique_matches_map) {
            for (const auto& bgc_pair : nrp_pair.second) {
                unique_matches.push_back(bgc_pair.second);
            }
        }

        // Step 1: Sort all matches by descending score.
        vector<MatchInfoLight> sortedMatches = std::move(unique_matches);
        sort(sortedMatches.begin(), sortedMatches.end(),
             [](const MatchInfoLight& a, const MatchInfoLight& b) {
            return a.score > b.score;
        });

        // Step 2: Filter by BGC. For each BGC (ignoring variant_idx), only keep up to max_num_matches_per_bgc.
        // Define BGC_key as a tuple of (genome_id, contig_idx, bgc_idx).
        unordered_map<BGC_Info, int> bgcCounts;
        vector<MatchInfoLight> bgcFiltered;
        for (const auto& match : sortedMatches) {
            auto bgcKey = match.bgc_info;
            variant_idx(bgcKey) = 0;  // Ignore variant_idx
            if (bgcCounts[bgcKey] < config.max_num_matches_per_bgc or config.max_num_matches_per_bgc == 0) {
                bgcFiltered.push_back(match);
                bgcCounts[bgcKey]++;
            }
        }

        // Step 3: Filter by NRP_ID. Only keep up to max_num_matches_per_nrp for each NRP.
        unordered_map<NRP_ID, int> nrpCounts;
        vector<MatchInfoLight> nrpFiltered;
        for (const auto& match : bgcFiltered) {
            if (nrpCounts[match.nrp_id] < config.max_num_matches_per_nrp or config.max_num_matches_per_nrp == 0) {
                nrpFiltered.push_back(match);
                nrpCounts[match.nrp_id]++;
            }
        }

        // Step 4: From the remaining matches, return only the top max_num_matches overall.
        if (nrpFiltered.size() > static_cast<size_t>(config.max_num_matches) and config.max_num_matches > 0) {
            nrpFiltered.resize(config.max_num_matches);
        }

        return nrpFiltered;
    }
}

MatchInfo get_full_match_info(const MatchInfoLight& match_light,
                              const HMM& hmm)
{
    MatchInfo fullMatch;
    fullMatch.score = 0;
    fullMatch.bgc_info = match_light.bgc_info;
    fullMatch.nrp_id = match_light.nrp_id;

    // For each linearization in the light match, compute its optimal path.
    for (const auto* linPtr : match_light.linearizations) {
        // Store a copy of the linearization.
        fullMatch.linearizations.push_back(*linPtr);
        // Compute the optimal HMM path using the monomer sequence from the linearization.
        // Here we assume that the initial state is 0 and the final state is the last state.
        auto [pathScore, optimalPath] = get_opt_hmm_path(hmm, 0, static_cast<StateIdx>(hmm.transitions.size() - 1),
                                                         linPtr->first);
        fullMatch.optimal_paths.push_back(optimalPath);
        fullMatch.score += pathScore;
    }

    return fullMatch;
}

vector<MatchInfo> get_matches(const unordered_map<BGC_Info, HMM>& hmms,
                              const vector<pair<NRP_Linearizations, NRP_ID>>& nrp_linearizations,
                              const MatchingConfig& config,
                              int num_threads)
{
    // Set the number of threads for OpenMP.
    omp_set_num_threads(num_threads);

    // Convert the unordered_map into a vector for easy index-based iteration.
    vector<pair<BGC_Info, HMM>> bgc_hmm_pairs(hmms.begin(), hmms.end());

    // This vector will accumulate MatchInfoLight objects from each BGC variant.
    vector<MatchInfoLight> matches_light;

    // Parallelize processing of each BGC variant using OpenMP.
#pragma omp parallel
    {
        // Each thread collects its local matches.
        vector<MatchInfoLight> local_matches;
#pragma omp for nowait
        for (int i = 0; i < static_cast<int>(bgc_hmm_pairs.size()); i++) {
            const auto& [bgc_info, hmm] = bgc_hmm_pairs[i];
            auto matches_light_for_bgc_variant = get_best_matches_for_bgc_variant(hmm, bgc_info,
                                                                                  nrp_linearizations,
                                                                                  config.max_num_matches_per_bgc);
            local_matches.insert(local_matches.end(),
                                 matches_light_for_bgc_variant.begin(),
                                 matches_light_for_bgc_variant.end());
        }
        // Merge the thread-local matches into the global vector.
#pragma omp critical
        {
            matches_light.insert(matches_light.end(), local_matches.begin(), local_matches.end());
        }
    }

    // Filter and sort the matches according to configuration limits.
    vector<MatchInfoLight> matches_light_filtered = filter_and_sort_matches(matches_light, config);

    // Reconstruct the full match info for each filtered match.
    vector<MatchInfo> matches;
    for (const auto& match_light : matches_light_filtered) {
        matches.push_back(get_full_match_info(match_light, hmms.at(match_light.bgc_info)));
    }
    return matches;
}