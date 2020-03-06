//
// Created by Maxwell Murphy on 3/4/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_COMBINATIONINDICESGENERATOR_H
#define TRANSMISSION_NETWORKS_APP_COMBINATIONINDICESGENERATOR_H

#include <vector>

// Adapted from https://stackoverflow.com/a/9432150/2755374

struct CombinationIndicesGenerator {
    using combination_t = std::vector<int>;

    bool completed;
    unsigned long generated = 1;

    CombinationIndicesGenerator(int n, int r);
    combination_t next() noexcept;

private:

    int n_;
    int r_;
    combination_t curr;
};


#endif //TRANSMISSION_NETWORKS_APP_COMBINATIONINDICESGENERATOR_H
