#include "parse_hmms.h"
#include <fstream>
#include <stdexcept>
#include "json.hpp" // https://github.com/nlohmann/json
#include <iostream>

std::unordered_map<BGC_Variant_ID , HMM>
parse_hmms_from_json(const std::string& hmm_json_path)
{
    using json = nlohmann::json;
    std::unordered_map<BGC_Variant_ID, HMM> hmms_map;

    std::ifstream ifs(hmm_json_path);
    if (!ifs.is_open()) {
        throw std::runtime_error("Failed to open HMM JSON file: " + hmm_json_path);
    }

    json j;
    ifs >> j;
    if (!j.is_array()) {
        throw std::runtime_error("HMM JSON top-level must be an array.");
    }

    for (auto& entry : j) {
        // Parse BGC_Info

        auto genome_id   = entry["bgc_variant_id"]["bgc_id"]["genome_id"].get<std::string>();
        int contig_idx   = entry["bgc_variant_id"]["bgc_id"]["contig_idx"].get<int>();
        int bgc_idx      = entry["bgc_variant_id"]["bgc_id"]["bgc_idx"].get<int>();
        int variant_idx  = entry["bgc_variant_id"]["variant_idx"].get<int>();
        BGC_ID bgc_id = std::make_tuple(genome_id, contig_idx, bgc_idx);
        BGC_Variant_ID bgc_variant_id = std::make_tuple(bgc_id, variant_idx);

        // Parse HMM
        HMM hmm;
        // transitions
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
        for (auto& em_row : entry["emissions"]) {
            std::vector<LogProb> em_vec;
            for (auto& val : em_row) {
                em_vec.push_back(val.get<LogProb>());
            }
            hmm.emissions.push_back(em_vec);
        }
        // module_start_states
        for (auto& val : entry["module_start_states"]) {
            hmm.module_start_states.push_back(val.get<StateIdx>());
        }
        // module_match_states
        for (auto& val : entry["module_match_states"]) {
            hmm.module_match_states.push_back(val.get<StateIdx>());
        }

        // Insert into map
        hmms_map[bgc_variant_id] = std::move(hmm);
    }

    return hmms_map;
}
