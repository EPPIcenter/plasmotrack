//
// Created by Maxwell Murphy on 4/24/23.
//

#include "gtest/gtest.h"

#include "core/parameters/Parameter.h"

#include "core/distributions/pdfs/DiscretePDF.h"


using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::distributions;
using namespace transmission_nets::core::computation;


TEST(DiscreteTest, CoreTest) {
    auto idx = std::make_shared<Parameter<int>>(0);
    std::vector<Probability> probs = {.1, .2, .3, .4};

    auto dist = std::make_shared<DiscreteDistribution>(probs);
    DiscretePDF pdf(idx, dist);

    ASSERT_DOUBLE_EQ(pdf.value(), std::log(.1));

    ASSERT_FALSE(pdf.isDirty());

    idx->saveState(1);
    idx->setValue(1);
    ASSERT_TRUE(pdf.isDirty());
    ASSERT_DOUBLE_EQ(pdf.value(), std::log(.2));
    idx->acceptState();
    ASSERT_FALSE(pdf.isDirty());

}