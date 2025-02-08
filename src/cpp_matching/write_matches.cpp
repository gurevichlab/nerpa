#include "write_matches.h"
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
                {"bgc_variant_info", {
                                  {"genome_id", std::get<0>(match.bgc_info)},
                                  {"contig_idx", std::get<1>(match.bgc_info)},
                                  {"bgc_idx", std::get<2>(match.bgc_info)},
                                  {"variant_idx", std::get<3>(match.bgc_info)}
                          }},
                {"nrp_id", match.nrp_id},
                {"nrp_linearizations", linearizations_json},
                {"optimal_paths", optimal_paths_json}
        };
    });

    std::ofstream ofs(output_path);
    if (!ofs.is_open()) {
        throw std::runtime_error("Failed to open output file: " + output_path);
    }
    ofs << j_out.dump(2) << std::endl;
}
