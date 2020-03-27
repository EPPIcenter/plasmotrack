//
// Created by Maxwell Murphy on 3/24/20.
//

#include "gtest/gtest.h"

#include "core/parameters/Parameter.h"

#include "model/transmission_process/node_transmission_process/GeometricGenerationProbability.h"


constexpr int MAX_TRANSMISSIONS = 20;

TEST(GeometricGenerationProbabilityTest, CoreTest) {
    Parameter<double> pGeoGenProb(.5);

    GeometricGenerationProbability<MAX_TRANSMISSIONS> geoGenProb(pGeoGenProb);

    ASSERT_DOUBLE_EQ(geoGenProb.value().sum(), 1);

    ASSERT_FALSE(geoGenProb.isDirty());
    pGeoGenProb.saveState();
    pGeoGenProb.setValue(.75);
    ASSERT_TRUE(geoGenProb.isDirty());
    ASSERT_DOUBLE_EQ(geoGenProb.value().sum(), 1);
    pGeoGenProb.acceptState();
    ASSERT_FALSE(geoGenProb.isDirty());

}