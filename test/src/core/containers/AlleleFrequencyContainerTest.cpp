//
// Created by Maxwell Murphy on 3/6/20.
//

#include "gtest/gtest.h"

#include "core/containers/AlleleFrequencyContainer.h"
#include "core/datatypes/Simplex.h"
#include "core/io/serialize.h"


#include <fmt/core.h>
#include <memory>
#include <vector>


using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::containers;
using namespace transmission_nets::core::parameters;

TEST(AlleleFrequencyContainerTest, HandlesChangedFrequencies) {
    using AlleleFrequenciesVector = Simplex;
    using AFC                     = AlleleFrequencyContainer<AlleleFrequenciesVector>;


    auto as1              = std::make_shared<Locus>("AS1", 4);
    auto as2              = std::make_shared<Locus>("AS2", 8);
    bool frequencyChanged = false;
    auto afc              = std::make_shared<AFC>();
    afc->addLocus(as1);
    afc->addLocus(as2);

    afc->alleleFrequencies(as1)->initializeValue({.4, .3, .2, .1});
    afc->alleleFrequencies(as2)->initializeValue({1, 1, 1, 1, 1, 1, 1, 1});

    fmt::print("Freq: {}\n", transmission_nets::core::io::serialize(afc->alleleFrequencies(as1)->value()));
    fmt::print("Freq: {}\n", transmission_nets::core::io::serialize(afc->alleleFrequencies(as2)->value()));

    afc->add_post_change_listener([&]() { frequencyChanged = true; });

    afc->alleleFrequencies(as1)->saveState(1);
    afc->alleleFrequencies(as1)->setValue({.1, .2, .3, .4});
    EXPECT_TRUE(frequencyChanged);
    afc->alleleFrequencies(as1)->restoreState(1);
}