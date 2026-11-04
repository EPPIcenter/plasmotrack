//
// Created by Maxwell Murphy on 12/9/19.
//

#include "core/datatypes/Genotype.h"
#include "core/parameters/Parameter.h"
#include "gtest/gtest.h"



using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::parameters;


class GenotypeTest : public ::testing::Test {
protected:
    GenotypeTest() : a1("111000111000"), a2("111111000000"), a3("10010"){};
    Genotype a1;
    Genotype a2;
    Genotype a3;
};

TEST_F(GenotypeTest, HandlesCounts) {
    ASSERT_EQ(a1.totalAlleles(), 12);
    ASSERT_EQ(Genotype::falsePositiveCount(a1, a2), 3);
    ASSERT_EQ(Genotype::falseNegativeCount(a1, a2), 3);
    ASSERT_EQ(Genotype::truePositiveCount(a1, a2), 3);
    ASSERT_EQ(Genotype::trueNegativeCount(a1, a2), 3);
    ASSERT_EQ(a1.totalPositiveCount(), 6);
    ASSERT_EQ(a1.totalAlleles(), 12);
    ASSERT_EQ(a3.totalAlleles(), 5);
}

TEST_F(GenotypeTest, HandlesFlipSetReset) {
    a1.flip(11);
    ASSERT_EQ(a1.totalPositiveCount(), 7);
    a1.flip(11);
    ASSERT_EQ(a1.totalPositiveCount(), 6);

    a1.set(8);
    ASSERT_EQ(a1.totalPositiveCount(), 6);
    a1.reset(8);
    ASSERT_EQ(a1.totalPositiveCount(), 5);
}

TEST_F(GenotypeTest, HandlesIteration) {
    // "10010"
    for (unsigned int i = 0; i < a3.totalAlleles(); ++i) {
        if (i == 0 or i == 3) {
            ASSERT_TRUE(a3.allele(i));
        } else {
            ASSERT_FALSE(a3.allele(i));
        }
    }
}

TEST_F(GenotypeTest, ParameterTest) {
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

TEST_F(GenotypeTest, InvertTest) {
    auto inverted_a1 = Genotype("000111000111");
    ASSERT_EQ(Genotype::invert(a1), inverted_a1);
    ASSERT_EQ(Genotype::invert(inverted_a1), a1);
    ASSERT_EQ(Genotype::invert(Genotype::invert(a1)), a1);
}

TEST_F(GenotypeTest, DiffTest) {
    Genotype a3_diff("10110");
    ASSERT_EQ(Genotype::diff(a3_diff, a3), Genotype("00100"));
}

TEST_F(GenotypeTest, SharedTest) {
    Genotype a3_shared("10100");
    ASSERT_EQ(Genotype::shared(a3_shared, a3), Genotype("10000"));
}

TEST_F(GenotypeTest, CopyTest) {
    Genotype a3_copy = a3;
    ASSERT_EQ(a3_copy, a3);
    a3.flip(1);
    ASSERT_NE(a3_copy, a3);
}

TEST_F(GenotypeTest, HandlesMutationMask) {
    Genotype g1("1011");
    Genotype g2("0011");
    Genotype g3("1010");
    Genotype g4("1001");

    ASSERT_EQ(g1.mutationMask(g2), Genotype("1000"));
    ASSERT_EQ(g1.mutationMask(g3), Genotype("0001"));
    ASSERT_EQ(g1.mutationMask(g4), Genotype("0010"));
    ASSERT_EQ(g1.mutationMask(g2).mutationMask(g3).mutationMask(g4), Genotype("0000"));
    ASSERT_EQ(g1.mutationMask(g2).mutationMask(g3), g1.mutationMask(g3).mutationMask(g2));
}

TEST_F(GenotypeTest, AValueBuiltFromAnAlleleCountInvertsOverItsAlleles) {
    Genotype last(4);// "0001"
    ASSERT_EQ(last, Genotype("0001"));
    EXPECT_EQ(Genotype::invert(last), Genotype("1110"));
}

TEST_F(GenotypeTest, SettingEveryAlleleStaysWithinTheAlleleCount) {
    Genotype g("0000");
    g.set();
    EXPECT_EQ(g, Genotype("1111"));
    EXPECT_EQ(g.totalPositiveCount(), 4u);
    EXPECT_EQ(Genotype::trueNegativeCount(g, g), 0u);
    EXPECT_EQ(g.allelesStr(), "1111");
}
