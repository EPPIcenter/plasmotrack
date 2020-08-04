//
// Created by Maxwell Murphy on 7/14/20.
//

//#include <cmath>
//
//#include "CombinationIndicesGenerator.h"
//#include "ProbAnyMissing.h"

/**
 * Calculate the probability that one or more events never occur over a sequence of trials
 * @param eventProbs Vector of probabilities of event A_1 ... A_n
 * @param numEvents Number of trials
 * @return
 */
//double probAnyMissing(const std::vector<double> eventProbs, int numEvents) {
//    int totalEvents = eventProbs.size();
//
//    if (numEvents < totalEvents) {
//        return 1.0;
//    }
//
//
//    double prob = 0;
//    double eventCombo;
//    CombinationIndicesGenerator c;
//    for (int i = 1; i <= totalEvents; ++i) {
//        c.reset(totalEvents, i);
//        while(!c.completed) {
//            eventCombo = 0.0;
//            for (const auto j : c.curr) {
//                eventCombo += eventProbs.at(j);
//            }
//            prob += pow((1 - eventCombo), numEvents) * pow(-1, (i - 1) % 2);
//            c.next();
//        }
//    }
//
//    return prob;
//}
