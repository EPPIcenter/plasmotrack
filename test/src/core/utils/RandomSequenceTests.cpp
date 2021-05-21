//
// Created by Maxwell Murphy on 4/17/20.
//

#include <boost/random.hpp>

#include "gtest/gtest.h"

#include "core/utils/generators/RandomSequence.h"

using namespace transmission_nets::core::utils::generators;

TEST(RandomSequenceTest, HandlesZeroIndexed) {

    boost::random::mt19937 r;
    auto test = randomSequence(0, 10, &r);
    ASSERT_EQ(test.size(), 10);
}
