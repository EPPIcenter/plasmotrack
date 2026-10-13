//
// Created by Maxwell Murphy on 7/14/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_PROBANYMISSING_H
#define TRANSMISSION_NETWORKS_APP_PROBANYMISSING_H

#include "core/computation/PartialLikelihood.h"
#include "core/utils/generators/CombinationIndicesGenerator.h"

#include <cmath>
#include <cstdint>
#include <vector>

// Architecture-specific intrinsics (x86/x64 only)
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
//#include <emmintrin.h>  // For SSE2 instructions
#include <immintrin.h>// For AVX instructions (if available)
#endif

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

        // The probability for each number of trials from 1 to numEvents. The result is overwritten by the next call.
        const std::vector<Likelihood>& vectorized(const std::vector<Likelihood>& eventProbs, unsigned int numEvents);


        // __m256 prob8{};
        // __m256 baseVec8{};
        // __m256 r8{};
        std::vector<Likelihood> baseVec{};
        std::vector<Likelihood> probVec{};
        // By n: the nonempty subsets of n events as bit masks, by size and then in lexicographic order.
        std::vector<std::vector<std::uint32_t>> subsetOrders{};
        const std::vector<std::uint32_t>& subsetOrder(std::size_t n);
        static std::size_t binomial(std::size_t n, std::size_t k);
        // std::vector<double> signVec{};
        generators::CombinationIndicesGenerator c;
    };

}// namespace transmission_nets::core::utils


#endif//TRANSMISSION_NETWORKS_APP_PROBANYMISSING_H
