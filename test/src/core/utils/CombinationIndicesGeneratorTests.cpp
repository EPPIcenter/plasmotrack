//
// Created by Maxwell Murphy on 3/6/20.
//

#include "core/utils/generators/CombinationIndicesGenerator.h"
#include "gtest/gtest.h"


#include <fmt/core.h>

using namespace transmission_nets::core::utils;

TEST(CombinationsIndicesGeneratorTest, HandlesGeneratingCombos) {
    generators::CombinationIndicesGenerator cs(10, 2);
    while (!cs.completed) {
        ASSERT_EQ(cs.curr.size(), 2);
        ASSERT_TRUE(cs.curr[0] < cs.curr[1]);
        cs.next();
    }
    ASSERT_TRUE(cs.completed);
    ASSERT_EQ(cs.generated, 45);
}

TEST(CombinationsIndicesGeneratorTest, HandlesGeneratingCombos2) {
    generators::CombinationIndicesGenerator cs(3, 2);
    while (!cs.completed) {
        ASSERT_EQ(cs.curr.size(), 2);
        ASSERT_TRUE(cs.curr[0] < cs.curr[1]);
        fmt::print("{} {}\n", cs.curr[0], cs.curr[1]);
        cs.next();
    }
    ASSERT_TRUE(cs.completed);
}

TEST(CombinationsIndicesGeneratorTest, HandlesGeneratingCombos3) {
    generators::CombinationIndicesGenerator cs(4, 1);
    while (!cs.completed) {
        ASSERT_EQ(cs.curr.size(), 1);
        cs.next();
    }
    ASSERT_TRUE(cs.completed);
}

TEST(CombinationsIndicesGeneratorTest, HandlesZeroChoice) {
    generators::CombinationIndicesGenerator cs(10, 0);
    ASSERT_EQ(cs.curr.size(), 0);
    ASSERT_TRUE(cs.completed);
}

TEST(CombinationsIndicesGeneratorTest, HandlesNoChoices) {
    generators::CombinationIndicesGenerator cs(0, 10);
    ASSERT_EQ(cs.curr.size(), 10);
    ASSERT_TRUE(cs.completed);
}