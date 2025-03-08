#include "write_matches.h"
#include "data_types.h"
#include <fstream>
#include <stdexcept>
#include "parsing/json.hpp"
#include <algorithm>

void write_matches_to_json(const std::vector<MatchInfo>& matches,
                           const std::string& output_path)
{
    using json = nlohmann::json;
    json j_out = json::array();

    std::transform(matches.begin(), matches.end(), std::back_inserter(j_out), [](const auto& match) {
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

    std::ofstream ofs(output_path);
    if (!ofs.is_open()) {
        throw std::runtime_error("Failed to open output file: " + output_path);
    }
    ofs << j_out.dump(2) << std::endl;
}
