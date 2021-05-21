//
// Created by Maxwell Murphy on 7/14/20.
//

#include "gtest/gtest.h"

#include "core/utils/ProbAnyMissing.h"

using namespace transmission_nets::core::utils;

TEST(ProbAnyMissingTests, TestSimpleVec) {

    probAnyMissingFunctor probAnyMissing;

    ASSERT_NEAR(probAnyMissing({.1, .2, .3, .4}, 5), 0.856, 1e-3);
    ASSERT_NEAR(probAnyMissing({.1, .2, .7}, 4), .832, 1e-3);

}
