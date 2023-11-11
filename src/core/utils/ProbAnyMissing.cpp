//
// Created by Maxwell Murphy on 7/14/20.
//

#include <vector>

#include "ProbAnyMissing.h"

namespace transmission_nets::core::utils {
    /**
     * Calculate the probability that one or more events never occur over a sequence of trials
     * @param eventProbs Vector of probabilities of event A_1 ... A_n
     * @param numEvents Number of trials
     * @return
     */

    Likelihood transmission_nets::core::utils::probAnyMissingFunctor::operator()(const std::vector<Likelihood>& eventProbs, unsigned int numEvents) {
                std::size_t totalEvents = eventProbs.size();
                int multCounter;
                Likelihood r;

                if (numEvents < totalEvents) {
                    return 1.0;
                }

                prob = 0.0;

                //      Calculate via inclusion-exclusion principle
                int sign = -1;

                baseVec.resize(0);
                for (std::size_t i = 1; i <= totalEvents; ++i) {
                    sign = -sign;

                    c.reset(totalEvents, i);

                    while (!c.completed) {
                        base = 1.0;
                        for (const auto j : c.curr) {
                            base -= eventProbs[j];
                        }

                        baseVec.push_back(base);
                        c.next();
                    }

                    // SIMD vectorization
                    if (baseVec.size() >= 4) {
                        for (std::size_t j = 0; j <= baseVec.size() - 4; j += 4) {
                            baseVec4 = _mm256_loadu_pd(&baseVec[j]);
                            r4       = _mm256_set_pd(sign, sign, sign, sign);

                            multCounter = (signed) numEvents;
                            while (multCounter > 0) {
                                if (multCounter & 1) {
                                    r4 = _mm256_mul_pd(r4, baseVec4);
                                }

                                baseVec4 = _mm256_mul_pd(baseVec4, baseVec4);
                                multCounter >>= 1;
                            }
                            prob += r4[0] + r4[1] + r4[2] + r4[3];
                        }
                    }

                    // scalar remainder
                    for (std::size_t j = baseVec.size() - (baseVec.size() % 4); j < baseVec.size(); ++j) {
//                        sign = -sign;
                        base = baseVec[j];
                        r    = sign;

                        multCounter = (signed) numEvents;
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


        //        std::size_t totalEvents = eventProbs.size();
        //        int multCounter;
        //        double r;
        //
        //        if (numEvents < totalEvents) {
        //            return 1.0;
        //        }
        //
        //        prob = 0.0;
        //
        //        //      Calculate via inclusion-exclusion principle
        //        int sign = -1;
        //        for (std::size_t i = 1; i <= totalEvents; ++i) {
        //            sign = -sign;
        //            c.reset(totalEvents, i);
        //            while (!c.completed) {
        //                base        = 1.0;
        //                multCounter = (signed) numEvents;
        //
        //                for (const auto j : c.curr) {
        //                    base -= eventProbs[j];
        //                }
        //                r = sign;
        //
        //                // squared exponentiation
        //                while (multCounter > 0) {
        //                    if (multCounter & 1) {
        //                        r *= base;
        //                    }
        //
        //                    base = (base * base);
        //                    multCounter >>= 1;
        //                }
        //
        //                prob += r;
        //                c.next();
        //            }
        //        }
        //
        //        return prob;
        //    }

//        std::size_t totalEvents = eventProbs.size();
//        int multCounter;
//        double r;
//
//        if (numEvents < totalEvents) {
//            return 1.0;
//        }
//
//        prob = 0.0;
//
//        //      Calculate via inclusion-exclusion principle
//        int sign = -1;
//        for (std::size_t i = 1; i <= totalEvents; ++i) {
//            sign = -sign;
//            c.reset(totalEvents, i);
//            baseVec.resize(0);
//            while (!c.completed) {
//                base = 1.0;
//
//                for (const auto j : c.curr) {
//                    base -= eventProbs[j];
//                }
//                baseVec.push_back(base);
//                c.next();
//            }
//
//            for (const double j : baseVec) {
//                base = j;
//
//                r = sign;
//
//                multCounter = (signed) numEvents;
//                // squared exponentiation
//                while (multCounter > 0) {
//                    if (multCounter & 1) {
//                        r *= base;
//                    }
//
//                    base = (base * base);
//                    multCounter >>= 1;
//                }
//
//                prob += r;
//            }
//        }
//
//        return prob;
//    }
}// namespace transmission_nets::core::utils
