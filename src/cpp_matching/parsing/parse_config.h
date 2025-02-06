#pragma once

#include <string>
#include "../data_types.h"

/**
 * Parses the matching configuration (max_num_matches_per_bgc, etc.) from JSON.
 * Example structure:
 * {
 *   "max_num_matches_per_bgc": 5,
 *   "max_num_matches_per_nrp": 10,
 *   "max_num_matches": 100
 * }
 *
 * @param config_json_path path to the JSON file
 * @return a MatchingConfig object
 */
MatchingConfig parse_config(const std::string& config_json_path);
