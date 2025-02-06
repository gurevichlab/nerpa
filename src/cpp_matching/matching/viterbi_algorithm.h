#ifndef NERPA2_VITERBI_ALGORITHM_H
#define NERPA2_VITERBI_ALGORITHM_H

#endif //NERPA2_VITERBI_ALGORITHM_H

#include "../data_types.h"

LogProb get_hmm_score(const HMM& hmm,
                      const vector<MonCode>& nrp_monomers,
                      const vector<pair<StateIdx, int>>& checkpoints);

pair<LogProb, vector<StateIdx>>
get_opt_hmm_path(const HMM& hmm,
                 StateIdx initial_state,
                 StateIdx final_state,
                 const vector<MonCode>& nrp_monomers);


