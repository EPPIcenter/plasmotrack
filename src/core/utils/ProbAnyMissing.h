//
// Created by Maxwell Murphy on 7/14/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_PROBANYMISSING_H
#define TRANSMISSION_NETWORKS_APP_PROBANYMISSING_H

#include "core/utils/generators/CombinationIndicesGenerator.h"
#include "core/computation/PartialLikelihood.h"

//#include <emmintrin.h>  // For SSE2 instructions
#include <immintrin.h>  // For AVX instructions (if available)

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
        Likelihood operator()(const std::vector<Likelihood>& eventProbs, unsigned int numEvents);

        __m256d prob4{};
        __m256d baseVec4{};
        __m256d r4{};
        double prob{};
        double base{};
        std::vector<double> baseVec{};
        
        generators::CombinationIndicesGenerator c;
    };

}// namespace transmission_nets::core::utils


#endif//TRANSMISSION_NETWORKS_APP_PROBANYMISSING_H
