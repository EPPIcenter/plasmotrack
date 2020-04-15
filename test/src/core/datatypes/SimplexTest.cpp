//
// Created by Maxwell Murphy on 3/9/20.
//

#include <vector>

#include "gtest/gtest.h"
#include "core/datatypes/Simplex.h"
#include "core/parameters/Parameter.h"

constexpr int MAX_ALLELES = 10;
TEST(AlleleFrequenciesVectorTest, BasicTest) {
    Simplex<MAX_ALLELES> av(3);
    ASSERT_DOUBLE_EQ(av.frequencies(0), 1);
    ASSERT_DOUBLE_EQ(av.alleleFrequencies().sum(), 1);

    av.set({.33, .33, .33});
    ASSERT_DOUBLE_EQ(av.frequencies(0), 1.0 / 3.0);
    ASSERT_DOUBLE_EQ(av.alleleFrequencies().sum(), 1);

    av.set({.25, .63, .33});
    ASSERT_DOUBLE_EQ(av.alleleFrequencies().sum(), 1);

    Simplex<MAX_ALLELES> av2({1, 2, 3});
    ASSERT_DOUBLE_EQ(av2.frequencies(0), 1.0 / 6.0);
    ASSERT_DOUBLE_EQ(av2.alleleFrequencies().sum(), 1);
}

TEST(AlleleFrequenciesVectorTest, ParameterTest) {
    Simplex<MAX_ALLELES> av(3);
    Parameter<Simplex<MAX_ALLELES>> p(av);
    ASSERT_DOUBLE_EQ(p.value().alleleFrequencies().sum(), 1);

    Parameter<Simplex<MAX_ALLELES>> p2({.33, .33, .33});
    ASSERT_DOUBLE_EQ(p2.value().alleleFrequencies(0), 1.0/3.0);
    ASSERT_DOUBLE_EQ(p2.value().alleleFrequencies().sum(), 1);

    p2.saveState();
    p2.setValue({.5, .6, .7});
    ASSERT_DOUBLE_EQ(p2.value().alleleFrequencies(0), .5 / (.5 + .6 + .7));
    p2.restoreState();
    ASSERT_DOUBLE_EQ(p2.value().alleleFrequencies(0), 1.0/3.0);

    p2.saveState();
    p2.setValue(Simplex<MAX_ALLELES>({1.0, 2.0, 3.0}));
    ASSERT_DOUBLE_EQ(p2.value().alleleFrequencies(0), 1.0 / 6.0);
    p2.restoreState();
    ASSERT_DOUBLE_EQ(p2.value().alleleFrequencies(0), 1.0/3.0);

    p2.saveState();
    p2.setValue(Simplex<MAX_ALLELES>(std::vector<double>({.1, .4, .6})));
    ASSERT_DOUBLE_EQ(p2.value().alleleFrequencies(0), .1 / (.1 + .4 + .6));
    p2.restoreState();
    ASSERT_DOUBLE_EQ(p2.value().alleleFrequencies(0), 1.0/3.0);
}