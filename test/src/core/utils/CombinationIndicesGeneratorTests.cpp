//
// Created by Maxwell Murphy on 3/6/20.
//

#include "core/utils/generators/CombinationIndicesGenerator.h"
#include "gtest/gtest.h"

using namespace transmission_nets::core::utils;

TEST(CombinationsIndicesGeneratorTest, HandlesGeneratingCombos) {
    generators::CombinationIndicesGenerator cs(10, 2);
    while(!cs.completed) {
        for (int i : cs.curr) {
            std::cout << i << ", ";
        }
        std::cout << std::endl;
        cs.next();
    }
    ASSERT_TRUE(cs.completed);
    ASSERT_EQ(cs.generated, 45);
}

TEST(CombinationsIndicesGeneratorTest, HandlesGeneratingCombos2) {
    generators::CombinationIndicesGenerator cs(4, 1);
    while(!cs.completed) {
        for (int i : cs.curr) {
            std::cout << i << ", ";
        }
        std::cout << std::endl;
        cs.next();
    }
    ASSERT_TRUE(cs.completed);
}

TEST(CombinationsIndicesGeneratorTest, HandlesZeroChoice) {
    generators::CombinationIndicesGenerator cs(10, 0);
    ASSERT_TRUE(cs.completed);
}

TEST(CombinationsIndicesGeneratorTest, HandlesNoChoices) {
   generators::CombinationIndicesGenerator cs(0, 10);
    ASSERT_TRUE(cs.completed);
}