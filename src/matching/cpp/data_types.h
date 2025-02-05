#include <vector>
using namespace std;


using MonCode = int;
using rBAN_idx = int;
using NRP_Linearization = pair<vector<MonCode>, vector<rBAN_idx>>;

class NRP_Linearizations {
public:
    vector<NRP_Linearization> non_iterative;
    vector<vector<vector<NRP_Linearization>>> iterative;
    // fragments are split into groups which are then linearized (and aligned) separately
    // linearizations of one group -> List[Linearization]
    // linearizations of all group members -> List[List[Linearization]]
    // linearizations of all groups -> List[List[List[Linearization]]]
};



using LogProb = float;
using StateIdx = int;

class HMM {
public:
    vector<vector<pair<StateIdx , LogProb>>> transitions;
    vector<vector<LogProb>> emissions;
    vector<vector<StateIdx>> nearest_module_start_state;
    // states are numbered from 0 to len(transitions) - 1
    // emissions[i][j] is the log probability of emitting MonCode j from state i
    // if state is not emitting, emissions[i] = {}
};

LogProb get_hmm_score(const HMM& hmm,
                      const vector<MonCode>& nrp_monomers,
                      const vector<pair<StateIdx, int>>& checkpoints);
