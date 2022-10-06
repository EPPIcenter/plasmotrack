//
// Created by Maxwell Murphy on 7/14/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_PROBANYMISSING_H
#define TRANSMISSION_NETWORKS_APP_PROBANYMISSING_H

#include "core/utils/generators/CombinationIndicesGenerator.h"
#include "core/computation/PartialLikelihood.h"

namespace transmission_nets::core::utils {

    using core::computation::Likelihood;
    struct probAnyMissingFunctor {

        probAnyMissingFunctor() = default;

        /**
         * Calculate the probability that one or more events never occur over a sequence of trials
         * @param eventProbs Vector of probabilities of event A_1 ... A_n
         * @param numEvents Number of trials
         * @return
         */
        double operator()(const std::vector<Likelihood>& eventProbs, unsigned int numEvents);

        Likelihood prob{};
        Likelihood eventCombo{};
        generators::CombinationIndicesGenerator c;
    };

}// namespace transmission_nets::core::utils


#endif//TRANSMISSION_NETWORKS_APP_PROBANYMISSING_H
