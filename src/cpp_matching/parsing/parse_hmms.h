#pragma once

#include <unordered_map>
#include <string>
#include "../data_types.h"


std::unordered_map<BGC_Variant_ID , HMM>
parse_hmms_from_json(const std::string& hmm_json_path);
