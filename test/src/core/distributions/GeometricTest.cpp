//
// Created by Maxwell Murphy on 3/24/20.
//

#include "gtest/gtest.h"

#include "core/parameters/Parameter.h"

#include "core/distributions/ZTGeometric.h"

using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::distributions;

constexpr int MAX_TRANSMISSIONS = 20;

TEST(GeometricTest, CoreTest) {
    auto pGeoGenProb = std::make_shared<Parameter<double>>(.5);
    //    Parameter<double> pGeoGenProb(.5);

    ZTGeometric<MAX_TRANSMISSIONS> geoProb(pGeoGenProb);

    ASSERT_DOUBLE_EQ(geoProb.value().sum(), 1);

    ASSERT_FALSE(geoProb.isDirty());
    pGeoGenProb->saveState(1);
    pGeoGenProb->setValue(.75);
    ASSERT_TRUE(geoProb.isDirty());
    ASSERT_DOUBLE_EQ(geoProb.value().sum(), 1);
    pGeoGenProb->acceptState();
    ASSERT_FALSE(geoProb.isDirty());
}