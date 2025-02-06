#include <iostream>
#include <string>
#include <stdexcept>
#include "data_types.h"
#include "parse_hmms.h"
#include "parse_nrps.h"
#include "parse_config.h"
#include "matcher.h"
#include "write_matches.h"

int main(int argc, char** argv)
{
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <hmms_json> <nrps_json> <config_json> <num_threads> [output_json]\n";
        return 1;
    }

    std::string hmms_json_path   = argv[1];
    std::string nrps_json_path   = argv[2];
    std::string config_json_path = argv[3];
    int num_threads              = std::stoi(argv[4]);
    std::string output_json_path = (argc > 5) ? argv[5] : "matches_out.json";

    try {
        // 1. Parse input data
        auto hmms_map           = parse_hmms_from_json(hmms_json_path);
        auto nrp_linearizations = parse_nrps_from_json(nrps_json_path);
        auto config             = parse_config(config_json_path);

        // 2. Perform matching
        auto matches = get_matches(hmms_map, nrp_linearizations, config, num_threads);

        // 3. Write results to JSON
        write_matches_to_json(matches, output_json_path);

        std::cout << "Done! Wrote " << matches.size() << " matches to " << output_json_path << "\n";
    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
    return 0;
}
