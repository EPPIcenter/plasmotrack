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
        for (std::size_t i = 1; i <= totalEvents; ++i) {
            sign = -sign;
            c.reset(totalEvents, i);
            while (!c.completed) {

                base = 1.0;
                multCounter = (signed) numEvents;

                for (const auto j : c.curr)
                {
                    base -= eventProbs[j];
                }


                r = sign;

                // squared exponentiation
                while (multCounter > 0)
                {
                    if (multCounter & 1)
                    {
                        r *= base;
                    }

                    base = (base * base);
                    multCounter >>= 1;
                }

                prob += r;
                c.next();
            }
        }

        return prob;
    }

}// namespace transmission_nets::core::utils
