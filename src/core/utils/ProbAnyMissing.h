//
// Created by Maxwell Murphy on 7/14/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_PROBANYMISSING_H
#define TRANSMISSION_NETWORKS_APP_PROBANYMISSING_H

#include <vector>
#include <cmath>

#include "CombinationIndicesGenerator.h"


//double probAnyMissing(const std::vector<double> eventProbs, int numEvents);


struct probAnyMissingFunctor {

    probAnyMissingFunctor() = default;

    /**
     * Calculate the probability that one or more events never occur over a sequence of trials
     * @param eventProbs Vector of probabilities of event A_1 ... A_n
     * @param numEvents Number of trials
     * @return
     */
    double operator()(const std::vector<double>& eventProbs, int numEvents) {
        int totalEvents = eventProbs.size();

        if (numEvents < totalEvents) {
            return 1.0;
        }

        prob = 0.0;

        for (int i = 1; i <= totalEvents; ++i) {
            c.reset(totalEvents, i);
            while(!c.completed) {
                eventCombo = 0.0;
                for (const auto j : c.curr) {
                    eventCombo += eventProbs.at(j);
                }
                prob += pow((1 - eventCombo), numEvents) * pow(-1, (i - 1) % 2);
                c.next();
            }
        }

        return prob;
    }

    double prob{};
    double eventCombo{};
    CombinationIndicesGenerator c;
};

#endif//TRANSMISSION_NETWORKS_APP_PROBANYMISSING_H
