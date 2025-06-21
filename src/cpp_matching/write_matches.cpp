#include "write_matches.h"
#include "data_types.h"
#include <fstream>
#include <stdexcept>
#include "parsing/json.hpp"
#include <algorithm>

void write_output_to_json(const std::vector<MatchInfo>& matches,
                          const std::unordered_map<BGC_Variant_ID, std::vector<double>>& p_values_vec,
                          const std::string& output_path)
{
    using json = nlohmann::json;

    // 1. Convert matches to JSON format
    json matches_json = json::array();

    std::transform(matches.begin(), matches.end(), std::back_inserter(matches_json), [](const auto& match) {
        json linearizations_json;
        std::transform(match.linearizations.begin(), match.linearizations.end(), std::back_inserter(linearizations_json),
                       [](const auto& lin) { return lin.second; });

        json optimal_paths_json;
        std::transform(match.optimal_paths.begin(), match.optimal_paths.end(), std::back_inserter(optimal_paths_json),
                       [](const auto& path) { return json(path); });


        return json{
                {"score", match.score},
                {"bgc_variant_id", {
                                  {"bgc_id", {
                                                     {"genome_id", genome_id(bgc_id(match.bgc_variant_id))},
                                                     {"contig_idx", contig_idx(bgc_id(match.bgc_variant_id))},
                                                     {"bgc_idx", bgc_idx(bgc_id(match.bgc_variant_id))}
                                             }},
                                  {"variant_idx", variant_idx(match.bgc_variant_id)}
                          }},
                {"nrp_id", match.nrp_id},
                {"nrp_linearizations", linearizations_json},
                {"optimal_paths", optimal_paths_json}
        }; // <-- Removed the extra semicolon here
    });

    // 2. Convert p-values to JSON format
    json p_values_json = json::array();
    for (const auto& [bgc_variant_id, pvals] : p_values_vec) {
        json bgc_variant_info = {
                {"bgc_id", {
                                   {"genome_id", genome_id(bgc_id(bgc_variant_id))},
                                   {"contig_idx", contig_idx(bgc_id(bgc_variant_id))},
                                   {"bgc_idx", bgc_idx(bgc_id(bgc_variant_id))}}},
                {"variant_idx", variant_idx(bgc_variant_id)}
        };
        p_values_json.push_back({
                                        {"bgc_variant_info", bgc_variant_info},
                                        {"p_values", pvals}
                                });
    }

    // 3. Combine matches and p-values into a single JSON object
    json j_out = {
            {"matches", matches_json},
            {"p_values", p_values_json}
    };

    // 4. Write the JSON object to the output file
    std::ofstream ofs(output_path);
    if (!ofs.is_open()) {
        throw std::runtime_error("Failed to open output file: " + output_path);
    }
    ofs << j_out.dump(2) << std::endl;
}
