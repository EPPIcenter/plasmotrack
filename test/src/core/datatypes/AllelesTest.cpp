//
// Created by Maxwell Murphy on 12/9/19.
//

#include "gtest/gtest.h"
#include "core/datatypes/Alleles.h"
#include "core/parameters/Parameter.h"

constexpr int MAX_ALLELES = 24;


using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::parameters;


using GeneticsImpl = AllelesBitSet<MAX_ALLELES>;


class AllelesTestFixture : public ::testing::Test {
protected:
    AllelesTestFixture() : a1("111000111000"), a2("111111000000"), a3("10010") {};
    GeneticsImpl a1;
    GeneticsImpl a2;
    GeneticsImpl a3;
};

TEST_F(AllelesTestFixture, HandlesCounts) {
    ASSERT_EQ(a1.totalAlleles(), 12);
    ASSERT_EQ(GeneticsImpl::falsePositiveCount(a1, a2), 3);
    ASSERT_EQ(GeneticsImpl::falseNegativeCount(a1, a2), 3);
    ASSERT_EQ(GeneticsImpl::truePositiveCount(a1, a2), 3);
    ASSERT_EQ(GeneticsImpl::trueNegativeCount(a1, a2), 3);
    ASSERT_EQ(a1.totalPositiveCount(), 6);
    ASSERT_EQ(a1.totalAlleles(), 12);
    ASSERT_EQ(a3.totalAlleles(), 5);
}

TEST_F(AllelesTestFixture, HandlesFlipSetReset) {
    a1.flip(11);
    ASSERT_EQ(a1.totalPositiveCount(), 7);
    a1.flip(11);
    ASSERT_EQ(a1.totalPositiveCount(), 6);

    a1.set(8);
    ASSERT_EQ(a1.totalPositiveCount(), 6);
    a1.reset(8);
    ASSERT_EQ(a1.totalPositiveCount(), 5);
}

TEST_F(AllelesTestFixture, HandlesIteration) {
    // "10010"
    for (unsigned int i = 0; i < a3.totalAlleles(); ++i) {
        if(i == 0 or i == 3) {
            ASSERT_TRUE(a3.allele(i));
        } else {
            ASSERT_FALSE(a3.allele(i));
        }
    }
}

TEST_F(AllelesTestFixture, ParameterTest) {
    Parameter<GeneticsImpl> p(a1);
    p.saveState("state1");
    p.setValue(GeneticsImpl("111000111001"));
    p.saveState("state2");
    p.setValue(GeneticsImpl("111000111011"));
    std::cout << p.value() << std::endl;
    p.restoreState("state2");
    std::cout << p.value() << std::endl;
    p.restoreState("state1");
    std::cout << p.value() << std::endl;
}
