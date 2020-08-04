//
// Created by Maxwell Murphy on 4/17/20.
//

#include <boost/random.hpp>

#include "core/utils/RandomSequence.h"
#include "gtest/gtest.h"


TEST(RandomSequenceTest, HandlesZeroIndexed) {

    boost::random::mt19937 r;
    auto test = randomSequence(0, 10, &r);
    ASSERT_EQ(test.size(), 10);
}
