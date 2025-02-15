#include <iostream>
#include <string>
#include <chrono>
#include "CLI11.hpp"

// Include the headers for your modules:
#include "parsing/parse_hmms.h"
#include "parsing/parse_nrps.h"
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
    int max_num_matches_per_bgc, max_num_matches_per_nrp, max_num_matches;
    std::string output_json_path;

    // Create command line options (very similar to argparse in Python)
    app.add_option("-H,--hmms_json", hmms_json_path,
                   "Path to the HMM/BGC JSON file")
            ->required();

    app.add_option("-N,--nrps_json", nrps_json_path,
                   "Path to the NRP linearizations JSON file")
            ->required();

    app.add_option("--max_num_matches_per_bgc", max_num_matches_per_bgc,
                   "Maximum number of matches to keep per BGC. Set to 0 to keep all matches.")
            ->required();

    app.add_option("--max_num_matches_per_nrp", max_num_matches_per_nrp,
                   "Maximum number of matches to keep per NRP. Set to 0 to keep all matches.")
            ->required();

    app.add_option("--max_num_matches", max_num_matches,
                   "Maximum number of matches to keep in total. Set to 0 to keep all matches.")
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
        std::cout << "Parsing NRP linearizations from " << hmms_json_path << std::endl;
        auto nrp_linearizations = parse_nrps_from_json(nrps_json_path);

        // 3. Parse the matching configuration
        auto matching_cfg = MatchingConfig{
            max_num_matches_per_bgc,
            max_num_matches_per_nrp,
            max_num_matches
        };

        // 4. Call your matching function to get the final list of matches
        std::cout << "Matching " <<
            hmms_map.size() << " HMMs against " <<
            nrp_linearizations.size() << " NRP linearizations in " <<
            num_threads << " threads...\n";

        auto start = std::chrono::high_resolution_clock::now(); // Start time
        auto matches = get_matches(hmms_map, nrp_linearizations, matching_cfg, num_threads);
        auto end = std::chrono::high_resolution_clock::now();   // End time
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Matching took " << duration.count() / 1000.0 << " s" << std::endl;

        // 5. Dump the resulting matches to an output JSON file
        std::cout << "Writing matches to " << output_json_path << std::endl;
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
