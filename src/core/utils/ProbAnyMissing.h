//
// Created by Maxwell Murphy on 7/14/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_PROBANYMISSING_H
#define TRANSMISSION_NETWORKS_APP_PROBANYMISSING_H

#include "CombinationIndicesGenerator.h"

namespace transmission_nets::core::utils {

    struct probAnyMissingFunctor {

        probAnyMissingFunctor() = default;

        /**
         * Calculate the probability that one or more events never occur over a sequence of trials
         * @param eventProbs Vector of probabilities of event A_1 ... A_n
         * @param numEvents Number of trials
         * @return
         */
        double operator()(const std::vector<double>& eventProbs, int numEvents);

        double prob{};
        double eventCombo{};
        CombinationIndicesGenerator c;
    };

}



#endif//TRANSMISSION_NETWORKS_APP_PROBANYMISSING_H
