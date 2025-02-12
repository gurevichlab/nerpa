#include "parse_config.h"
#include <fstream>
#include <stdexcept>
#include "json.hpp"

MatchingConfig parse_config(const std::string& config_json_path)
{
    using json = nlohmann::json;
    std::ifstream ifs(config_json_path);
    if (!ifs.is_open()) {
        throw std::runtime_error("Failed to open config JSON file: " + config_json_path);
    }

    json j;
    ifs >> j;

    MatchingConfig config;
    config.max_num_matches_per_bgc = j.value("max_num_matches_per_bgc", 0);
    config.max_num_matches_per_nrp = j.value("max_num_matches_per_nrp", 0);
    config.max_num_matches         = j.value("max_num_matches", 0);

    return config;
}
