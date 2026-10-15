//
// Created by Maxwell Murphy on 12/9/19.
//

#include "core/datatypes/Genotype.h"
#include "core/parameters/Parameter.h"
#include "gtest/gtest.h"

constexpr int MAX_ALLELES = 24;


using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::parameters;


class AllelesTestFixture : public ::testing::Test {
protected:
    AllelesTestFixture() : a1("111000111000"), a2("111111000000"), a3("10010"){};
    Genotype a1;
    Genotype a2;
    Genotype a3;
};

TEST_F(AllelesTestFixture, HandlesCounts) {
    ASSERT_EQ(a1.totalAlleles(), 12);
    ASSERT_EQ(Genotype::falsePositiveCount(a1, a2), 3);
    ASSERT_EQ(Genotype::falseNegativeCount(a1, a2), 3);
    ASSERT_EQ(Genotype::truePositiveCount(a1, a2), 3);
    ASSERT_EQ(Genotype::trueNegativeCount(a1, a2), 3);
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
        if (i == 0 or i == 3) {
            ASSERT_TRUE(a3.allele(i));
        } else {
            ASSERT_FALSE(a3.allele(i));
        }
    }
}

TEST_F(AllelesTestFixture, ParameterTest) {
    Parameter<Genotype> p(a1);
    p.saveState(1);
    p.setValue(Genotype("111000111001"));
    p.saveState(2);
    p.setValue(Genotype("111000111011"));
    std::cout << p.value() << std::endl;
    p.restoreState(2);
    std::cout << p.value() << std::endl;
    p.restoreState(1);
    std::cout << p.value() << std::endl;
}

TEST_F(AllelesTestFixture, InvertTest) {
    auto inverted_a1 = Genotype("000111000111");
    ASSERT_EQ(Genotype::invert(a1), inverted_a1);
    ASSERT_EQ(Genotype::invert(inverted_a1), a1);
    ASSERT_EQ(Genotype::invert(Genotype::invert(a1)), a1);
}

TEST_F(AllelesTestFixture, DiffTest) {
    Genotype a3_diff("10110");
    ASSERT_EQ(Genotype::diff(a3_diff, a3), Genotype("00100"));
}

TEST_F(AllelesTestFixture, SharedTest) {
    Genotype a3_shared("10100");
    ASSERT_EQ(Genotype::shared(a3_shared, a3), Genotype("10000"));
}

TEST_F(AllelesTestFixture, CopyTest) {
    Genotype a3_copy = a3;
    ASSERT_EQ(a3_copy, a3);
    a3.flip(1);
    ASSERT_NE(a3_copy, a3);
}

TEST(AllelesTest, HandlesMutationMask) {
    Genotype a1("1011");
    Genotype a2("0011");
    Genotype a3("1010");
    Genotype a4("1001");

    ASSERT_EQ(a1.mutationMask(a2), Genotype("1000"));
    ASSERT_EQ(a1.mutationMask(a3), Genotype("0001"));
    ASSERT_EQ(a1.mutationMask(a4), Genotype("0010"));
    ASSERT_EQ(a1.mutationMask(a2).mutationMask(a3).mutationMask(a4), Genotype("0000"));
    ASSERT_EQ(a1.mutationMask(a2).mutationMask(a3), a1.mutationMask(a3).mutationMask(a2));


}
