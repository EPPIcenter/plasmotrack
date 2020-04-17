//
// Created by Maxwell Murphy on 3/6/20.
//

#include "gtest/gtest.h"

#include "core/datatypes/Simplex.h"
#include "core/containers/AlleleFrequencyContainer.h"


TEST(AlleleFrequencyContainerTest, HandlesChangedFrequencies) {
    using AlleleFrequenciesVector = Simplex;
    using AlleleFrequencyContainer = AlleleFrequencyContainer<AlleleFrequenciesVector>;

    Locus as1("AS1", 4);
    Locus as2("AS2", 8);

    bool frequencyChanged = false;

    std::vector<AlleleFrequencyContainer::LocusAlleleFrequencyAssignment> lfas {
            {&as1, AlleleFrequenciesVector{.4, .3, .2, .1}},
            {&as2, AlleleFrequenciesVector{1, 1, 1, 1, 1, 1, 1, 1}}
    };

    AlleleFrequencyContainer afc(lfas);
    afc.add_post_change_listener([&]() { frequencyChanged = true; });

    std::cout << afc << std::endl;
    afc.alleleFrequencies(&as1).saveState();
    afc.alleleFrequencies(&as1).setValue({.1, .2, .3, .4});
    EXPECT_TRUE(frequencyChanged);
    afc.alleleFrequencies(&as1).restoreState();
}