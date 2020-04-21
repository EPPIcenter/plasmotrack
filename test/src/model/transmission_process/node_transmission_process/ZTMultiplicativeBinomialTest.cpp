//
// Created by Maxwell Murphy on 4/20/20.
//

#include "gtest/gtest.h"

#include "core/parameters/Parameter.h"

#include "model/transmission_process/node_transmission_process/ZTMultiplicativeBinomial.h"

constexpr int MAX_COI = 15;

TEST(ZTMultiplicativeBinomialTest, CoreTest) {
    Parameter<double> prob(.5);
    Parameter<double> assoc(1.0);

    ZTMultiplicativeBinomial<MAX_COI> ztmb(prob, assoc);

    ASSERT_TRUE(ztmb.isDirty());
    ASSERT_FALSE(ztmb.value().hasNaN());
    ASSERT_FALSE(ztmb.isDirty());
}