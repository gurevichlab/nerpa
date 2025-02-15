#pragma once
#include "../data_types.h"
#include <vector>

LogProb get_hmm_score(const HMM& hmm,
                      const std::vector<MonCode>& nrp_monomers,
                      const std::vector<std::pair<StateIdx, int>>& checkpoints);

std::pair<LogProb, std::vector<StateIdx>>
get_opt_hmm_path(const HMM& hmm,
                 StateIdx initial_state,
                 StateIdx final_state,
                 const std::vector<MonCode>& nrp_monomers);


