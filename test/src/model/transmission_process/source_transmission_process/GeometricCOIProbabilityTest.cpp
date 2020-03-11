//
// Created by Maxwell Murphy on 3/10/20.
//

#include "gtest/gtest.h"
#include "model/transmission_process/source_transmission_process/GeometricCOIProbability.h"

constexpr int MAX_COI = 10;
TEST(GeometricCOIProbabilityTest, BasicTest) {
    Parameter<double> prob(.3);
    GeometricCOIProbability<MAX_COI> coiProb(prob);
    EXPECT_DOUBLE_EQ(coiProb.value().sum(), 1);

    auto tmp = coiProb.value();

    prob.saveState();
    prob.setValue(.6);
    EXPECT_TRUE(coiProb.isDirty());
    EXPECT_DOUBLE_EQ(coiProb.value().sum(), 1);
    EXPECT_FALSE(coiProb.isDirty());
    prob.restoreState();

    for (int i = 0; i < tmp.size(); ++i) {
        EXPECT_EQ(tmp(i), coiProb.value()(i));
    }
}
