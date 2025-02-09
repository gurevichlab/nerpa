#include "parse_nrps.h"
#include <fstream>
#include <stdexcept>
#include "json.hpp"

std::vector<std::pair<NRP_Linearizations, NRP_ID>>
parse_nrps_from_json(const std::string& nrp_json_path)
{
    using json = nlohmann::json;
    std::vector<std::pair<NRP_Linearizations, NRP_ID>> nrp_linearizations;

    std::ifstream ifs(nrp_json_path);
    if (!ifs.is_open()) {
        throw std::runtime_error("Failed to open NRP JSON file: " + nrp_json_path);
    }

    json j;
    ifs >> j;
    if (!j.is_array()) {
        throw std::runtime_error("NRP JSON top-level must be an array.");
    }

    for (auto& entry : j) {
        NRP_ID nrp_id = entry["nrp_id"].get<std::string>();

        NRP_Linearizations nrp_obj;
        // Parse non_iterative
        for (auto& lin_j : entry["non_iterative"]) {
            NRP_Linearization lin;
            auto mon_codes_j = lin_j[0];
            auto rban_idxs_j = lin_j[1];

            for (auto& m : mon_codes_j) {
                lin.first.push_back(m.get<MonCode>());
            }
            for (auto& r : rban_idxs_j) {
                lin.second.push_back(r.get<rBAN_idx>());
            }
            nrp_obj.non_iterative.push_back(lin);
        }

        // Parse iterative
        for (auto& split_j : entry["iterative"]) {
            std::vector<std::vector<NRP_Linearization>> split_vec;
            for (auto& group_j : split_j) {
                std::vector<NRP_Linearization> group_vec;
                for (auto& lin_j : group_j) {
                    NRP_Linearization lin;
                    auto mon_codes_j = lin_j[0];
                    auto rban_idxs_j = lin_j[1];
                    for (auto& m : mon_codes_j) {
                        lin.first.push_back(m.get<MonCode>());
                    }
                    for (auto& r : rban_idxs_j) {
                        lin.second.push_back(r.get<rBAN_idx>());
                    }
                    group_vec.push_back(lin);
                }
                split_vec.push_back(group_vec);
            }
            nrp_obj.iterative.push_back(split_vec);
        }

        nrp_linearizations.emplace_back(std::move(nrp_obj), nrp_id);
    }

    return nrp_linearizations;
}
