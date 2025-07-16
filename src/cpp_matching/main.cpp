#include <iostream>
#include <string>
#include <chrono>
#include "CLI11.hpp"
#include "openmp_wrapper.h"

// Include the headers for your modules:
#include "parsing/parse_hmms.h"
#include "parsing/parse_nrps.h"
#include "matching/matcher.h"
#include "p_values/p_values.h"
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
        auto parsing_results = parse_hmms_from_json(hmms_json_path);
        auto &hmms_for_matching = parsing_results.first;
        auto &hmms_for_p_value_estimation = parsing_results.second;
        //auto [hmms_for_matching, hmms_for_p_value_estimation] = parse_hmms_from_json(hmms_json_path);

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
                  hmms_for_matching.size() << " HMMs against " <<
                  nrp_linearizations.size() << " NRP linearizations in " <<
            num_threads << " threads...\n";

        auto start = std::chrono::high_resolution_clock::now(); // Start time
        auto matches = get_matches(hmms_for_matching, nrp_linearizations, matching_cfg, num_threads);
        auto end = std::chrono::high_resolution_clock::now();   // End time
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Matching took " << duration.count() / 1000.0 << " s" << std::endl;

        // 5. Estimate p-values for the HMMs
        std::cout << "Estimating p-values in " << num_threads << " threads..." << std::endl;
        start = std::chrono::high_resolution_clock::now(); // Start time

        std::unordered_map<BGC_Variant_ID, std::vector<double>> p_values_by_hmm;
        p_values_by_hmm.reserve(hmms_for_matching.size());  // Reserve space for efficiency
        std::vector<BGC_Variant_ID> bgc_variant_ids;
        bgc_variant_ids.reserve(hmms_for_matching.size());
        for (const auto& key_value : hmms_for_matching) {
            bgc_variant_ids.push_back(key_value.first);
        }

#pragma omp parallel
        {
            std::vector<std::pair<BGC_Variant_ID, std::vector<double>>> p_values_for_thread;
#pragma omp for nowait
            for (int i = 0; i < static_cast<int>(bgc_variant_ids.size()); ++i) {
                const auto &bgc_variant_id = bgc_variant_ids[i];
                const auto &hmm_lo = hmms_for_matching[bgc_variant_id];
                const auto &hmm_lp = hmms_for_p_value_estimation.at(bgc_variant_id);
                auto p_values = compute_hmm_p_values(hmm_lo, hmm_lp);
                p_values_for_thread.emplace_back(bgc_variant_id, p_values);
            }
#pragma omp critical
            {
                for (const auto &[bgc_variant_id, p_values]: p_values_for_thread) {
                    p_values_by_hmm[bgc_variant_id] = std::move(p_values);
                }
            }
        }

        end = std::chrono::high_resolution_clock::now();   // End time
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Estimated p-values for " << p_values_by_hmm.size() << " HMMs." << std::endl;
        std::cout << "P-value estimation took " << duration.count() / 1000.0 << " s" << std::endl;


        // 6. Dump the resulting matches to an output JSON file
        std::cout << "Writing matches and p-values to " << output_json_path << std::endl;
        write_output_to_json(matches, p_values_by_hmm, output_json_path);

        std::cout << "Done! Wrote " << matches.size() << " matches to "
                  << output_json_path << "\n";
    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
