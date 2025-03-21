#pragma once

#include <string>
#include <vector>
#include "data_types.h"

/**
 * Writes the final vector of MatchInfo to a JSON file.
 *
 * @param matches The vector of MatchInfo to write.
 * @param output_path The path to the JSON file to create or overwrite.
 */
void write_matches_to_json(const std::vector<MatchInfo>& matches,
                           const std::string& output_path);
