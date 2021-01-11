//
// Created by Maxwell Murphy on 3/6/20.
//

#include "gtest/gtest.h"
#include "core/utils/CombinationIndicesGenerator.h"

using namespace transmission_nets::core::utils;

TEST(CombinationsIndicesGeneratorTest, HandlesGeneratingCombos) {
    CombinationIndicesGenerator cs(10, 2);
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
    CombinationIndicesGenerator cs(4, 1);
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
    CombinationIndicesGenerator cs(10, 0);
    ASSERT_TRUE(cs.completed);
}

TEST(CombinationsIndicesGeneratorTest, HandlesNoChoices) {
    CombinationIndicesGenerator cs(0, 10);
    ASSERT_TRUE(cs.completed);
}