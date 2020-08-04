//
// Created by Maxwell Murphy on 7/14/20.
//

#include "gtest/gtest.h"

#include "core/utils/ProbAnyMissing.h"


TEST(ProbAnyMissingTests, TestSimpleVec) {

    std::vector<double> pvec = {.1, .2, .3, .4};
    int n_trials = 5;
    probAnyMissingFunctor probAnyMissing;

    double prob = probAnyMissing(pvec, n_trials);
    ASSERT_NEAR(prob, 0.856, 1e-3);

}
