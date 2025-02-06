#include "write_matches.h"
#include <fstream>
#include <stdexcept>
#include <nlohmann/json.hpp>

void write_matches_to_json(const std::vector<MatchInfo>& matches,
                           const std::string& output_path)
{
    using json = nlohmann::json;
    json j_out = json::array();

    for (const auto& match : matches) {
        json match_j;
        match_j["score"] = match.score;

        // BGC_Info is a tuple<string,int,int,int>
        match_j["bgc_info"] = {
                std::get<0>(match.bgc_info),
                std::get<1>(match.bgc_info),
                std::get<2>(match.bgc_info),
                std::get<3>(match.bgc_info)
        };

        match_j["nrp_id"] = match.nrp_id;

        // Convert linearizations
        {
            json lin_array = json::array();
            for (const auto& lin : match.linearizations) {
                json lin_j;
                lin_j["first"]  = lin.first;
                lin_j["second"] = lin.second;
                lin_array.push_back(lin_j);
            }
            match_j["linearizations"] = lin_array;
        }

        // Convert optimal_paths
        {
            json paths_array = json::array();
            for (const auto& path : match.optimal_paths) {
                paths_array.push_back(path);
            }
            match_j["optimal_paths"] = paths_array;
        }

        j_out.push_back(match_j);
    }

    std::ofstream ofs(output_path);
    if (!ofs.is_open()) {
        throw std::runtime_error("Failed to open output file: " + output_path);
    }
    ofs << j_out.dump(2) << std::endl;
}
