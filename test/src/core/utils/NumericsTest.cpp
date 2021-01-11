//
// Created by Maxwell Murphy on 9/21/20.
//

#include "gtest/gtest.h"

#include "core/utils/numerics.h"

using namespace transmission_nets::core::utils;

TEST(NumericsTests, TestExpNormalize) {
    std::vector<double> input{1, -10, 1000};
    auto out = expNormalize(input);
    std::vector<double> expected{0, 0, 1};

    for (size_t i = 0; i < out.size(); ++i) {
        std::cout << "Input/Output - " << input[i] << " " << out[i] << std::endl;
        ASSERT_NEAR(out[i], expected[i], 1e-6);
    }


}