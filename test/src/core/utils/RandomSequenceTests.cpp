//
// Created by Maxwell Murphy on 4/17/20.
//


#include "gtest/gtest.h"

#include "core/io/serialize.h"
#include "core/utils/generators/RandomSequence.h"

#include <boost/random.hpp>
#include <fmt/core.h>

#include <algorithm>

using namespace transmission_nets::core::utils::generators;
using namespace transmission_nets::core::io;

TEST(RandomSequenceTest, HandlesZeroIndexed) {
    auto r    = std::make_shared<boost::random::mt19937>();
    auto test = randomSequence(0, 10, r);
    fmt::print("Random Sequence: {}\n", serialize(test));
    ASSERT_EQ(*std::min_element(test.begin(), test.end()), 0);
    ASSERT_EQ(*std::max_element(test.begin(), test.end()), 9);
    ASSERT_EQ(test.size(), 10);
}
