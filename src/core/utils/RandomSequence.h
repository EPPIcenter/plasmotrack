//
// Created by Maxwell Murphy on 4/17/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_RANDOMSEQUENCE_H
#define TRANSMISSION_NETWORKS_APP_RANDOMSEQUENCE_H

#include <vector>
#include <numeric>

#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <boost/range/algorithm.hpp>

/**
 * Generate a random sequence of values from [min, max). Useful for randomly indexing into a vector
 * @tparam Engine boost random generator
 * @param min  min inclusive value in sequence
 * @param max max exclusive value in sequence
 * @param rng boost random generator
 * @return vector<int> of random sequence containing [min, max)
 */
template<typename Engine>
std::vector<int> randomSequence(int min, int max, Engine rng) {
    assert(min < max);
    std::vector<int> indices(max - min);
    std::iota(std::begin(indices), std::end(indices), min);

    auto int_generator = [=](int max_val) {
        boost::random::uniform_int_distribution<> dist_{0, max_val};
        return dist_(*rng);
    };

    boost::range::random_shuffle(indices, int_generator);
    return indices;
}

#endif//TRANSMISSION_NETWORKS_APP_RANDOMSEQUENCE_H
