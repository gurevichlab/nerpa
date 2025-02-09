#include <iostream>
#include <string>
#include "CLI11.hpp"

// Include the headers for your modules:
#include "parsing/parse_hmms.h"
#include "parsing/parse_nrps.h"
#include "parsing/parse_config.h"
#include "matching/matcher.h"
#include "write_matches.h"

int main(int argc, char** argv)
{
    CLI::App app{"NRP Matching Tool"};

    // Define variables to store the parsed arguments
    std::string hmms_json_path;
    std::string nrps_json_path;
    std::string config_json_path;
    int num_threads;
    std::string output_json_path;

    // Create command line options (very similar to argparse in Python)
    app.add_option("-H,--hmms_json", hmms_json_path,
                   "Path to the HMM/BGC JSON file")
            ->required();

    app.add_option("-N,--nrps_json", nrps_json_path,
                   "Path to the NRP linearizations JSON file")
            ->required();

    app.add_option("-C,--config_json", config_json_path,
                   "Path to the matching configuration JSON file")
            ->required();

    app.add_option("-t,--threads", num_threads,
                   "Number of OpenMP threads to use (e.g., 4)")
            ->required();

    app.add_option("-o,--output", output_json_path,
                   "Output JSON file for the resulting matches")
            ->default_val("matches_out.json");

    // Let CLI11 parse the command line
    CLI11_PARSE(app, argc, argv);

    try {
        // 1. Parse input data (HMMs, BGC_Info)
        std::cout << "Parsing HMMs from " << hmms_json_path << std::endl;
        auto hmms_map = parse_hmms_from_json(hmms_json_path);

        // 2. Parse NRP linearizations + NRP IDs
        auto nrp_linearizations = parse_nrps_from_json(nrps_json_path);

        // 3. Parse the matching configuration
        auto config = parse_config(config_json_path);

        // 4. Call your matching function to get the final list of matches
        auto matches = get_matches(hmms_map, nrp_linearizations, config, num_threads);

        // 5. Dump the resulting matches to an output JSON file
        write_matches_to_json(matches, output_json_path);

        std::cout << "Done! Wrote " << matches.size() << " matches to "
                  << output_json_path << "\n";
    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
