#pragma once

#include <unordered_map>
#include <vector>
#include "../data_types.h"

std::vector<MatchInfo> get_matches(const std::unordered_map<BGC_Variant_ID, HMM>& hmms,
                                   const std::vector<std::pair<NRP_Linearizations, NRP_ID>>& nrp_linearizations,
                                   const MatchingConfig& config,
                                   int num_threads);
