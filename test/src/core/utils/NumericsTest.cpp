//
// Created by Maxwell Murphy on 9/21/20.
//

#include "gtest/gtest.h"

#include "core/utils/numerics.h"

#include <fmt/core.h>

using namespace transmission_nets::core::utils;

TEST(NumericsTests, TestExpNormalize) {
    std::vector<double> input{1, -10, 1000};
    auto out = expNormalize(input);
    std::vector<double> expected{0, 0, 1};

    for (size_t i = 0; i < out.size(); ++i) {
        ASSERT_NEAR(out[i], expected[i], 1e-6);
    }
}

TEST(NumericsTests, TestLogitVector) {
    std::vector<double> input{.001, .01, .1};
    auto out = logit(input);
    std::vector<double> expected{-6.906755, -4.595120, -2.197225};

    for (size_t i = 0; i < out.size(); ++i) {
        ASSERT_NEAR(out[i], expected[i], 1e-6);
    }
}

TEST(NumericsTests, TestLogitSum) {
    std::vector<double> input{-1.31192, -3.66517, -1.83803, -1.94771, -1.38826, -2.51401, -4.25089, -3.66517, -2.5140};
    auto out        = logitSum(input);
    double expected = -0.1191479;
    ASSERT_NEAR(out, expected, 1e-6);
}

TEST(NumericsTests, TestLogitScale) {
    std::vector<double> input{-1.31192, -3.66517, -1.83803, -1.94771, -1.38826, -2.51401, -4.25089, -3.66517, -2.51401};
    double scale = -0.280057;
    auto out     = logitScale(input, scale);
    std::vector<double> expected{-1.655684, -3.951461, -2.156220, -2.262006, -1.727477, -2.813645, -4.534422, -3.951461, -2.813645};

    for (size_t i = 0; i < out.size(); ++i) {
        ASSERT_NEAR(out[i], expected[i], 1e-6);
    }
}

TEST(NumericsTests, TestInitKVec) {
    auto out = initKvecs<4, 3>();
    for (const auto & kvec : out) {
        for (const auto & k : kvec) {
           fmt::print("{} ", k);
        }
        fmt::print("\n");
    }

}