//
// Created by Maxwell Murphy on 4/20/20.
//

#include "gtest/gtest.h"

#include "core/parameters/Parameter.h"

#include "core/distributions/ZTMultiplicativeBinomial.h"

using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::distributions;

constexpr int MAX_COI = 15;

TEST(ZTMultiplicativeBinomialTest, CoreTest) {
    auto prob  = std::make_shared<Parameter<double>>(.5);
    auto assoc = std::make_shared<Parameter<double>>(1.0);

    ZTMultiplicativeBinomial<MAX_COI> ztmb(prob, assoc);

    ASSERT_FALSE(ztmb.isDirty());
    ASSERT_FALSE(ztmb.value().hasNaN());
    ASSERT_FALSE(ztmb.isDirty());
}