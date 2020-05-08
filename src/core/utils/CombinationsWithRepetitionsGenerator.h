//
// Created by Maxwell Murphy on 4/29/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_COMBINATIONSWITHREPETITIONSGENERATOR_H
#define TRANSMISSION_NETWORKS_APP_COMBINATIONSWITHREPETITIONSGENERATOR_H

#include <vector>
#include <numeric>

struct CombinationsWithRepetitionsGenerator {
    using combination_t = std::vector<int>;

    bool completed = false;
    unsigned long generated = 0;

    CombinationsWithRepetitionsGenerator(unsigned int nChoices, unsigned int k);
    void next() noexcept;
    combination_t curr{};

private:

    int nChoices_;
    int k_;
};


#endif //TRANSMISSION_NETWORKS_APP_COMBINATIONSWITHREPETITIONSGENERATOR_H
