#include "data_types.h"
#include <vector>
#include <limits>
#include <utility>
using namespace std;

// Assuming these types are defined somewhere:
using LogProb = double;

// Forward declarations for functions used in lambdas.
// double alignment_score(const AlignmentType& alignment);
// double get_hmm_score(const HMM& hmm, const std::vector<int>& seq);
LogProb linearization_score(const HMM& hmm,
                            const NRP_Linearization& linearization){
    vector<pair<StateIdx, int>> checkpoints = { {0, 0},
                                                {(int) hmm.transitions.size() - 1, linearization.first.size()} };
    return get_hmm_score(hmm, linearization.first, checkpoints);
}


std::pair<LogProb , std::vector<NRP_Linearization*>> get_best_linearizations_for_nrp(
        const HMM& hmm,
        const NRP_Linearizations& nrp_linearizations)
{

    // Determine the best non-iterative linearization.
    LogProb best_non_iterative_score = -std::numeric_limits<double>::infinity();
    NRP_Linearization* best_non_iterative_linearization;  // Default constructed; assumes non-empty nrp_linearizations.non_iterative

    for (const auto& non_iterative_linearization : nrp_linearizations.non_iterative) {
        LogProb score = linearization_score(non_iterative_linearization);
        if (score > best_non_iterative_score) {
            best_non_iterative_score = score;
            best_non_iterative_linearization = &non_iterative_linearization;
        }
    }

    // Determine the best iterative linearizations.
    double best_iterative_score = -std::numeric_limits<double>::infinity();
    std::vector<NRP_Linearization*> best_iterative_linearizations;

    // Assuming nrp_linearizations.iterative is a vector of "groups",
    // where each group is a vector (of groups) of linearizations.
    for (const auto& linearizations_of_one_split : nrp_linearizations.iterative) {
        double split_score = 0;
        std::vector<NRP_Linearization*> split_linearizations;
        for (const auto& group_linearizations : linearizations_of_one_split) {
            double group_score = -std::numeric_limits<double>::infinity();
            NRP_Linearization* best_group_linearization;
            for (const auto& linearization : group_linearizations) {
                double score = linearization_score(linearization);
                if (score > group_score) {
                    group_score = score;
                    best_group_linearization = &linearization;
                }
            }
            split_linearizations.push_back(best_group_linearization);
            split_score += group_score;
        }
        if (split_score > best_iterative_score) {
            best_iterative_score = split_score;
            best_iterative_linearizations = split_linearizations;
        }
    }


    // Return the best overall result.
    if (best_non_iterative_score > best_iterative_score) {
        return { best_non_iterative_score, { best_non_iterative_linearization } };
    } else {
        return { best_iterative_score, best_iterative_linearizations };
    }
}