//
// Created by Maxwell Murphy on 4/29/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_COMBINATIONSWITHREPETITIONSGENERATOR_H
#define TRANSMISSION_NETWORKS_APP_COMBINATIONSWITHREPETITIONSGENERATOR_H

#include <vector>

namespace transmission_nets::core::utils::generators {

    struct CombinationsWithRepetitionsGenerator {
        using combination_t = std::vector<int>;

        bool completed          = false;
        unsigned long generated = 0;

        CombinationsWithRepetitionsGenerator(int nChoices, int k);
        CombinationsWithRepetitionsGenerator();

        void next() noexcept;
        void reset(int nChoices, int k) noexcept;

        combination_t curr{};

    private:
        int nChoices_;
        int k_;
    };

}// namespace transmission_nets::core::utils::generators


#endif//TRANSMISSION_NETWORKS_APP_COMBINATIONSWITHREPETITIONSGENERATOR_H
