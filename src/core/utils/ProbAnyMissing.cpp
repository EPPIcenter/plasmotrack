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

    namespace {
        // Appends every r-subset of {first, ..., n - 1} added to `prefix`, as a bit mask, with subsets in lexicographic
        // order of their indices.
        void appendSubsets(const std::size_t n, const std::size_t first, const std::size_t r, const std::uint32_t prefix, std::vector<std::uint32_t>& out) {
            if (r == 0) {
                out.push_back(prefix);
                return;
            }
            for (std::size_t j = first; j + r <= n; ++j) {
                appendSubsets(n, j + 1, r - 1, prefix | (1u << j), out);
            }
        }
    }// namespace

    const std::vector<std::uint32_t>& probAnyMissingFunctor::subsetOrder(const std::size_t n) {
        if (subsetOrders.size() <= n) {
            subsetOrders.resize(n + 1);
        }
        auto& order = subsetOrders[n];
        if (order.empty()) {
            for (std::size_t i = 1; i <= n; ++i) {
                appendSubsets(n, 0, i, 0, order);
            }
        }
        return order;
    }

    const std::vector<Likelihood>& probAnyMissingFunctor::vectorized(const std::vector<Likelihood>& eventProbs, unsigned int numEvents) {
        const std::size_t totalEvents = eventProbs.size();

        probVec.assign(numEvents, 0.0);

        if (numEvents < totalEvents) {
            std::fill_n(probVec.begin(), numEvents, 1.0);
            return probVec;
        }

        std::fill_n(probVec.begin(), totalEvents - 1, 1.0);

        // 1 - p[c_0] - p[c_1] - ... for every subset, subtracting in index order: a subset's value is that of the
        // subset without its largest index, less that index's probability.
        baseVec.resize(std::size_t{1} << totalEvents);
        baseVec[0] = 1.0;
        for (std::size_t j = 0; j < totalEvents; ++j) {
            const std::size_t top = std::size_t{1} << j;
            for (std::size_t mask = 0; mask < top; ++mask) {
                baseVec[top | mask] = baseVec[mask] - eventProbs[j];
            }
        }

        //      Calculate via inclusion-exclusion principle, over subsets by size and then in lexicographic order
        const auto& order = subsetOrder(totalEvents);
        std::size_t next  = 0;
        double sign = -1.0;
        for (std::size_t i = 1; i <= totalEvents; ++i) {
            sign = -sign;
            const std::size_t count = binomial(totalEvents, i);
            for (std::size_t n = 0; n < count; ++n) {
                const Likelihood base = baseVec[order[next++]];
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

    std::size_t probAnyMissingFunctor::binomial(const std::size_t n, const std::size_t k) {
        std::size_t out = 1;
        for (std::size_t i = 1; i <= k; ++i) {
            out = out * (n - k + i) / i;
        }
        return out;
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
