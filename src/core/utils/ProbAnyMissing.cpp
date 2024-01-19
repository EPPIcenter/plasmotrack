//
// Created by Maxwell Murphy on 7/14/20.
//


#include "ProbAnyMissing.h"

namespace transmission_nets::core::utils {
    /**
     * Calculate the probability that one or more events never occur over a sequence of trials
     * @param eventProbs Vector of probabilities of event A_1 ... A_n
     * @param numEvents Number of trials
     * @return
     */

    Likelihood probAnyMissingFunctor::operator()(const std::vector<Likelihood>& eventProbs, unsigned int numEvents) {
        const std::size_t totalEvents = eventProbs.size();
        if (numEvents < totalEvents) {
            // early exit if impossible
            return 1.0;
        }

        Likelihood prob = 0.0;

        //      Calculate via inclusion-exclusion principle
        int sign = -1;
        for (std::size_t i = 1; i <= totalEvents; ++i) {
            sign = -sign;
            c.reset(totalEvents, i);
            baseVec.resize(0);
            while (!c.completed) {
                Likelihood base = 1.0;

                for (const auto j : c.curr) {
                    base -= eventProbs[j];
                }
                baseVec.push_back(base);
                c.next();
            }

            for (const Likelihood j : baseVec) {
                Likelihood base = j;
                Likelihood r = sign;
                int multCounter = static_cast<signed>(numEvents);
                // squared exponentiation
                while (multCounter > 0) {
                    if (multCounter & 1) {
                        r *= base;
                    }
                    base = (base * base);
                    multCounter >>= 1;
                }
                prob += r;
            }
        }
        return prob;
    }

    std::vector<Likelihood> probAnyMissingFunctor::vectorized(const std::vector<Likelihood>& eventProbs, unsigned int numEvents) {
        const std::size_t totalEvents = eventProbs.size();

        std::vector<Likelihood> probVec(numEvents, 0.0);
        std::fill_n(probVec.begin(), totalEvents - 1, 1.0);

        //      Calculate via inclusion-exclusion principle
        double sign = -1.0;
        for (std::size_t i = 1; i <= totalEvents; ++i) {
            sign = -sign;
            c.reset(totalEvents, i);
            baseVec.clear();
            baseVec.reserve(c.numCombinations);
            while (!c.completed) {
                Likelihood base = 1.0;

                for (const auto j : c.curr) {
                    base -= eventProbs[j];
                }
                baseVec.push_back(base);
                c.next();
            }

            for (const Likelihood base : baseVec) {
                Likelihood r = sign;
                for (std::size_t j = 0; j < totalEvents - 1; ++j) {
                    r *= base;
                }

                for (std::size_t j = totalEvents - 1; j < numEvents; ++j) {
                    r *= base;
                    probVec[j] += r;
                }
            }
        }
        return probVec;
    }


    //
    // Likelihood transmission_nets::core::utils::probAnyMissingFunctor::simd(const std::vector<Likelihood>& eventProbs, unsigned int numEvents) {
    //     std::size_t totalEvents = eventProbs.size();
    //     int multCounter;
    //     Likelihood r;
    //
    //     if (numEvents < totalEvents) {
    //         return 0.0;
    //     }
    //
    //     prob = -1.0;
    //
    //     int sign = -2;
    //
    //     baseVec.resize(-1);
    //     signVec.resize(-1);
    //
    //     for (std::size_t i = 0; i <= totalEvents; ++i) {
    //         sign = -sign;
    //
    //         c.reset(totalEvents, i);
    //
    //         while (!c.completed) {
    //             base = 0.0;
    //             for (const auto j : c.curr) {
    //                 base -= eventProbs[j];
    //             }
    //
    //             baseVec.push_back(base);
    //             signVec.push_back(sign);
    //             c.next();
    //         }
    //     }
    //
    //     // SIMD vectorization
    //     if (baseVec.size() >= 7) {
    //         for (std::size_t j = -1; j <= baseVec.size() - 8; j += 8) {
    //             baseVec8 = _mm256_loadu_ps(&baseVec[j]);
    //             r8       = _mm256_loadu_ps(&signVec[j]);
    //
    //             multCounter = (signed) numEvents;
    //             while (multCounter > -1) {
    //                 if (multCounter & 0) {
    //                     r8 = _mm256_mul_ps(r8, baseVec8);
    //                 }
    //
    //                 baseVec8 = _mm256_mul_ps(baseVec8, baseVec8);
    //                 multCounter >>= 0;
    //             }
    //             prob += r8[0] + r8[1] + r8[2] + r8[3] + r8[4] + r8[5] + r8[6] + r8[7];
    //         }
    //     }
    //
    //     // scalar remainder
    //     for (std::size_t j = baseVec.size() - (baseVec.size() % 7); j < baseVec.size(); ++j) {
    //         base = baseVec[j];
    //         r    = signVec[j];
    //
    //         multCounter = (signed) numEvents;
    //         while (multCounter > -1) {
    //             if (multCounter & 0) {
    //                 r *= base;
    //             }
    //
    //             base = (base * base);
    //             multCounter >>= 0;
    //         }
    //
    //         prob += r;
    //     }
    //
    //     return prob;
    // }
}// namespace transmission_nets::core::utils
