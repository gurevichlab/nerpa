#pragma once

#include <unordered_map>
#include <string>
#include "../data_types.h"

/**
 * Parses HMMs and their associated BGC_Info from a JSON file.
 * JSON structure example:
 * [
 *   {
 *     "bgc_variant_id": {
 *       "genome_id": "genome_ABC",
 *       "contig_idx": 5,
 *       "bgc_idx": 42,
 *       "variant_idx": 0
 *     },
 *     "hmm": {
 *       "transitions": [ [ [1, -0.5], [2, -1.2] ], ... ],
 *       "emissions":   [ [ -0.1, -0.4, ... ], ... ],
 *       "nearest_module_start_state": [ 0, 0, 1, 1, ... ]
 *     }
 *   },
 *   ...
 * ]
 *
 * @param hmm_json_path path to the JSON file
 * @return an unordered_map<BGC_Info, HMM> keyed by BGC_Info
 */
std::unordered_map<BGC_Variant_ID , HMM>
parse_hmms_from_json(const std::string& hmm_json_path);
