#pragma once

#include <vector>
#include <string>
#include <utility>
#include "data_types.h"

/**
 * Parses NRP linearizations from a JSON file.
 * Example JSON structure:
 * [
 *   {
 *     "nrp_id": "nrp_ABC",
 *     "non_iterative": [
 *       {
 *         "first": [1, 3, 2, 4],
 *         "second": [10, 11]
 *       },
 *       ...
 *     ],
 *     "iterative": [
 *       [  // first "split"
 *         [ // first "group" within this split
 *           { "first": [...], "second": [...] },
 *           ...
 *         ],
 *         [ // second "group" in this split
 *           ...
 *         ]
 *       ],
 *       ...
 *     ]
 *   },
 *   ...
 * ]
 *
 * @param nrp_json_path path to the JSON file
 * @return vector of (NRP_Linearizations, NRP_ID)
 */
std::vector<std::pair<NRP_Linearizations, NRP_ID>>
parse_nrps_from_json(const std::string& nrp_json_path);
