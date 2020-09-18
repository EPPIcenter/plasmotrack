//
// Created by Maxwell Murphy on 4/29/20.
//

#include "gtest/gtest.h"
#include "core/utils/CombinationsWithRepetitionsGenerator.h"

using namespace transmission_nets::core::utils;

TEST(CombinationsWithRepetitionsTest, CoreTest) {
    CombinationsWithRepetitionsGenerator cs(4, 3);

    std::vector<int> t{1,1,1,1,0,0};
    while (!cs.completed) {
        auto tmp = t;
//        auto indices = cs.next();
        cs.next();
        for (auto i : cs.curr) {
            tmp[i]++;
        }
//        for (auto j : tmp) {
//            std::cout << j << " ";
//        }
//        std::cout << "\n";
    }
//    std::cout << "Total Generated: " << cs.generated << std::endl;
    ASSERT_EQ(cs.generated, 20);
}

TEST(CombinationsWithRepetitionsTest, CoreTest2) {
    CombinationsWithRepetitionsGenerator cs(1, 3);
    while (!cs.completed) {
//        auto indices = cs.next();
        cs.next();
//        for (auto i : indices) {
//            std::cout << i << " ";
//        }
//        std::cout << "\n";
    }
//    std::cout << "Total Generated: " << cs.generated << std::endl;
    ASSERT_EQ(cs.generated, 1);
}

TEST(CombinationsWithRepetitionsTest, CoreTest3) {
    CombinationsWithRepetitionsGenerator cs(3, 1);
    while (!cs.completed) {
//        auto indices = cs.next();
        cs.next();
        for (auto i : cs.curr) {
            std::cout << i << " ";
        }
        std::cout << "\n";
    }
//    std::cout << "Total Generated: " << cs.generated << std::endl;
    ASSERT_EQ(cs.generated, 3);
}