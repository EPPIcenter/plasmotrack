//
// Created by Maxwell Murphy on 3/4/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_COMBINATIONINDICESGENERATOR_H
#define TRANSMISSION_NETWORKS_APP_COMBINATIONINDICESGENERATOR_H

#include <vector>


namespace transmission_nets::core::utils::generators {
    // Adapted from https://stackoverflow.com/a/9432150/2755374
    struct CombinationIndicesGenerator {
        using combination_t = std::vector<unsigned int>;

        bool completed;
        unsigned long generated = 1;
        unsigned long numCombinations = 0;
        combination_t curr{};

        /**
         * Generate a sequences of indices representing n choose r element combinations.
         * @param n number of elements
         * @param r number of choices
         */
        CombinationIndicesGenerator(std::size_t n, std::size_t r) noexcept;

        CombinationIndicesGenerator() noexcept;

        void reset(std::size_t n, std::size_t r) noexcept;
        void reset() noexcept;

        void next() noexcept;
        void advance(int n) noexcept;

    private:
        std::size_t n_;
        std::size_t r_;

        void calculateNumCombinations() noexcept;

    };
}// namespace transmission_nets::core::utils::generators


#endif//TRANSMISSION_NETWORKS_APP_COMBINATIONINDICESGENERATOR_H
