#include "parse_hmms.h"
#include <fstream>
#include <stdexcept>
#include <nlohmann/json.hpp> // https://github.com/nlohmann/json

std::unordered_map<BGC_Info, HMM>
parse_hmms_from_json(const std::string& hmm_json_path)
{
    using json = nlohmann::json;
    std::unordered_map<BGC_Info, HMM> hmms_map;

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
        auto genome_id   = entry["bgc_info"]["genome_id"].get<std::string>();
        int contig_idx   = entry["bgc_info"]["contig_idx"].get<int>();
        int bgc_idx      = entry["bgc_info"]["bgc_idx"].get<int>();
        int variant_idx  = entry["bgc_info"]["variant_idx"].get<int>();
        BGC_Info bgc_info = std::make_tuple(genome_id, contig_idx, bgc_idx, variant_idx);

        // Parse HMM
        HMM hmm;
        // transitions
        for (auto& st : entry["hmm"]["transitions"]) {
            std::vector<std::pair<StateIdx, LogProb>> row;
            for (auto& pair_j : st) {
                StateIdx st_to = pair_j[0].get<StateIdx>();
                LogProb  prob  = pair_j[1].get<LogProb>();
                row.emplace_back(st_to, prob);
            }
            hmm.transitions.push_back(row);
        }
        // emissions
        for (auto& em_row : entry["hmm"]["emissions"]) {
            std::vector<LogProb> em_vec;
            for (auto& val : em_row) {
                em_vec.push_back(val.get<LogProb>());
            }
            hmm.emissions.push_back(em_vec);
        }
        // nearest_module_start_state
        for (auto& val : entry["hmm"]["nearest_module_start_state"]) {
            hmm.nearest_module_start_state.push_back(val.get<StateIdx>());
        }

        // Insert into map
        hmms_map[bgc_info] = std::move(hmm);
    }

    return hmms_map;
}
