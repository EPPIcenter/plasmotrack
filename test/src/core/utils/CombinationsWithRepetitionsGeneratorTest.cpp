//
// Created by Maxwell Murphy on 4/29/20.
//

#include "core/utils/generators/CombinationsWithRepetitionsGenerator.h"
#include "gtest/gtest.h"

using namespace transmission_nets::core::utils::generators;

TEST(CombinationsWithRepetitionsTest, CoreTest) {
    CombinationsWithRepetitionsGenerator cs(4, 3);

    std::vector<int> t{1,1,1,1,0,0};
    while (!cs.completed) {
        auto tmp = t;
        cs.next();
        for (auto i : cs.curr) {
            tmp[i]++;
        }
   }
    ASSERT_EQ(cs.generated, 20);
}

TEST(CombinationsWithRepetitionsTest, CoreTest2) {
    CombinationsWithRepetitionsGenerator cs(1, 3);
    while (!cs.completed) {
        cs.next();
    }
    ASSERT_EQ(cs.generated, 1);
}

TEST(CombinationsWithRepetitionsTest, CoreTest3) {
    CombinationsWithRepetitionsGenerator cs(3, 1);
    while (!cs.completed) {
        cs.next();
        for (auto i : cs.curr) {
            std::cout << i << " ";
        }
        std::cout << "\n";
    }
    ASSERT_EQ(cs.generated, 3);
}