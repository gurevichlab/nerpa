#include "parse_hmms.h"
#include <fstream>
#include <stdexcept>
#include "json.hpp" // https://github.com/nlohmann/json
#include <iostream>
#include <limits>

HMM parse_hmm(const nlohmann::json& entry)
{
    HMM hmm;

    // parse BGC Variant ID ID
    // std::cout << "Parsing BGC Variant ID from JSON entry." << std::endl;
    auto genome_id   = entry["bgc_variant_id"]["bgc_id"]["genome_id"].get<std::string>();
    int contig_idx   = entry["bgc_variant_id"]["bgc_id"]["contig_idx"].get<int>();
    int bgc_idx      = entry["bgc_variant_id"]["bgc_id"]["bgc_idx"].get<int>();
    int variant_idx  = entry["bgc_variant_id"]["variant_idx"].get<int>();
    BGC_ID bgc_id = std::make_tuple(genome_id, contig_idx, bgc_idx);
    BGC_Variant_ID bgc_variant_id = std::make_tuple(bgc_id, variant_idx);

    hmm.bgc_variant_id = bgc_variant_id;

    // transitions
    // std::cout << "Parsing transitions from JSON entry." << std::endl;
    for (auto& st : entry["transitions"]) {
        std::vector<std::pair<StateIdx, LogProb>> row;
        for (auto& pair_j : st) {
            StateIdx st_to = pair_j[0].get<StateIdx>();
            LogProb  prob  = pair_j[1].get<LogProb>();
            row.emplace_back(st_to, prob);
        }
        hmm.transitions.push_back(row);
    }
    // emissions
    // std::cout << "Parsing emissions from JSON entry." << std::endl;
    for (auto& em_row : entry["emissions"]) {
        std::vector<LogProb> em_vec;
        for (auto& val : em_row) {
            if (val.is_null()) {
                // null → -∞
                em_vec.push_back(-std::numeric_limits<LogProb>::infinity());
            } else {
                // normal finite value
                em_vec.push_back(val.get<LogProb>());
            }
        }
        hmm.emissions.push_back(em_vec);
    }
    // module_start_states
    // std::cout << "Parsing module start states from JSON entry." << std::endl;
    for (auto& val : entry["module_start_states"]) {
        hmm.module_start_states.push_back(val.get<StateIdx>());
    }
    // module_match_states
    // std::cout << "Parsing module match states from JSON entry." << std::endl;
    for (auto& val : entry["module_match_states"]) {
        hmm.module_match_states.push_back(val.get<StateIdx>());
    }

    return hmm;
}

// 1. hmms for matching; 2. hmms for p-value estimation
std::pair<
        std::unordered_map<BGC_Variant_ID , HMM>,
        std::unordered_map<BGC_Variant_ID , HMM>
        >
parse_hmms_from_json(const std::string& hmm_json_path)
{
    using json = nlohmann::json;
    std::unordered_map<BGC_Variant_ID, HMM> hmms_for_matching;
    std::unordered_map<BGC_Variant_ID, HMM> hmms_for_p_value_estimation;

    std::ifstream ifs(hmm_json_path);
    if (!ifs.is_open()) {
        throw std::runtime_error("Failed to open HMM JSON file: " + hmm_json_path);
    }

    json j;
    ifs >> j;
    if (!j.is_array()) {
        throw std::runtime_error("HMM JSON top-level must be an array.");
    }

    for (auto& hmm_info : j) {
        // std::cout << "Parsing HMM for matching" << std::endl;
        auto hmm_for_matching = parse_hmm(hmm_info["hmm_for_matching"]);
        // std::cout << "Parsing HMM for p-value estimation" << std::endl;
        auto hmm_for_p_value_estimation = parse_hmm(hmm_info["hmm_for_p_values_estimation"]);
        auto bgc_variant_id = hmm_for_matching.bgc_variant_id;
        assert (bgc_variant_id == hmm_for_p_value_estimation.bgc_variant_id &&
               "BGC_Variant_IDs for matching and p-value estimation must match.");

        // Insert into maps
        hmms_for_matching[bgc_variant_id] = std::move(hmm_for_matching);
        hmms_for_p_value_estimation[bgc_variant_id] = std::move(hmm_for_p_value_estimation);
    }

    return std::make_pair(hmms_for_matching,
                          hmms_for_p_value_estimation);
}
