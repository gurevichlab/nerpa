#pragma once

#include <unordered_map>
#include <vector>
#include "../data_types.h"

/**
 * The main function that orchestrates matching.
 *
 * @param hmms An unordered_map of BGC_Info -> HMM.
 * @param nrp_linearizations A vector of (NRP_Linearizations, NRP_ID).
 * @param config The matching configuration.
 * @param num_threads Number of threads to use for parallelization.
 * @return A vector of MatchInfo with full HMM path reconstructions.
 */
std::vector<MatchInfo> get_matches(const std::unordered_map<BGC_Info, HMM>& hmms,
                                   const std::vector<std::pair<NRP_Linearizations, NRP_ID>>& nrp_linearizations,
                                   const MatchingConfig& config,
                                   int num_threads);
