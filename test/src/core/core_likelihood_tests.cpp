//
// Created by Maxwell Murphy on 12/9/19.
//

#include <iostream>

#include "gtest/gtest.h"
#include "core/Parameters/SimpleParameter.h"
#include "core/Likelihood/AdderLikelihood.h"
#include "core/Likelihood/Likelihood.h"
#include "core/Parameters/Alleles.h"

TEST(CoreLikelihoodTest, LikelihoodTest) {
    IntegerParameter a("a", 3);
    IntegerParameter b("b", 4);
    IntegerParameter c("c", 5);
    ASSERT_EQ(a.value(), 3);
    ASSERT_EQ(b.value(), 4);
    ASSERT_EQ(c.value(), 5);

    AdderLikelihood adder_like1(&a, &b);
    AdderLikelihood adder_like2(&a, &c);
    Likelihood lik;

    lik.pushLikelihood(&adder_like1);
    lik.pushLikelihood(&adder_like2);
    lik.reinitialize();
    ASSERT_EQ(lik.value(), 15);


    a.saveState();
    a.setValue(5);
    ASSERT_EQ(a.value(), 5);
    ASSERT_EQ(lik.value(), 19);
    a.acceptState();

    a.saveState();
    b.saveState();
    a.setValue(6);
    b.setValue(6);
    ASSERT_TRUE(adder_like1.isDirty());
    ASSERT_TRUE(adder_like2.isDirty());
    ASSERT_TRUE(lik.isDirty());
    ASSERT_EQ(a.value(), 6);
    ASSERT_EQ(b.value(), 6);
    ASSERT_EQ(lik.value(), 23);
    a.restoreState();
    b.restoreState();
    ASSERT_FALSE(adder_like1.isDirty());
    ASSERT_FALSE(adder_like2.isDirty());
    ASSERT_FALSE(lik.isDirty());
    ASSERT_FALSE(a.isDirty());
    ASSERT_FALSE(b.isDirty());
    ASSERT_FALSE(c.isDirty());
    ASSERT_EQ(a.value(), 5);
    ASSERT_EQ(b.value(), 4);
    ASSERT_EQ(c.value(), 5);
    ASSERT_EQ(lik.value(), 19);

    Alleles<32> test_alleles("test1", "011101");
    std::cout << test_alleles.id() << "\n";
    std::cout << test_alleles.value() << "\n";
    std::cout << test_alleles.num_alleles() << "\n";
}
