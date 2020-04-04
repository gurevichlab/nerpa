//
// Created by olga on 02.03.19.
//

#ifndef NRPSMATCHER_SCOREPRISM_H
#define NRPSMATCHER_SCOREPRISM_H

#include "Matcher/Score/Base/Score.h"

namespace matcher {
    class ScorePrism : public Score {
    public:
        bool getScore(const aminoacid::Aminoacid &nrpAA, const aminoacid::Aminoacid &predAA,
                        const nrpsprediction::AAdomainPrediction::AminoacidProb &prob,
                        const std::pair<int, int> &pos,
                        double& score) const override;
    };
}


#endif //NRPSMATCHER_SCOREPRISM_H
