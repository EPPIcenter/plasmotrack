//
// Created by mmurphy on 7/7/24.
//
// File: ZTGeometricTest.cpp

#include "gtest/gtest.h"
#include "core/distributions/ZTGeometric.h"

using namespace transmission_nets::core::distributions;

class ZTGeometricTest : public ::testing::Test {
protected:
    void SetUp() override {
        prob = std::make_shared<Parameter<double>>(0.5);
        ztGeometric = std::make_shared<ZTGeometric<10>>(prob);
    }

    std::shared_ptr<Parameter<double>> prob;
    std::shared_ptr<ZTGeometric<10>> ztGeometric;
};

TEST_F(ZTGeometricTest, InitializationTest) {
    ASSERT_NE(ztGeometric, nullptr);
}

TEST_F(ZTGeometricTest, ValueTest) {
    auto result = ztGeometric->value();

    double sum = 0.0;
    for (int i = 1; i <= 10; ++i) {
        sum += result(i);
    }

    ASSERT_NEAR(sum, 1.0, 1e-5) << "The sum of probabilities should be approximately 1.0";
}

TEST_F(ZTGeometricTest, ProbabilityUpdateTest) {
    prob->setValue(0.2);
    auto result = ztGeometric->value();

    double sum = 0.0;
    for (int i = 1; i <= 10; ++i) {
        sum += result(i);
    }

    ASSERT_NEAR(sum, 1.0, 1e-5) << "The sum of probabilities should be approximately 1.0 after updating probability";
}

// int main(int argc, char **argv) {
//     ::testing::InitGoogleTest(&argc, argv);
//     return RUN_ALL_TESTS();
// }