//
// Created by Maxwell Murphy on 3/5/20.
//

#include "gtest/gtest.h"
#include "core/parameters/Ordering.h"

TEST(OrderingTest, HandlesSwapsNotifies) {
    auto el1 = 1;
    auto el2 = 2;
    auto el3 = 3;
    auto el4 = 4;
    bool el1MovedRight = false;
    bool el2MovedLeft = false;

    Ordering<int> ordering;
    ordering.addElements({&el1, &el2, &el3, &el4});
    ordering.add_keyed_moved_left_listener(&el1, [&](const int* el) {
        el2MovedLeft = (el == &el2);
    });
    ordering.add_keyed_moved_right_listener(&el2, [&](const int* el) {
        el1MovedRight = (el == &el1);
    });

    ASSERT_EQ(ordering.value()[0], &el1);
    ASSERT_EQ(ordering.value()[1], &el2);
    ASSERT_EQ(el1, 1);
    ASSERT_EQ(el2, 2);
    ordering.swap(0, 1);

    ASSERT_EQ(ordering.value()[0], &el2);
    ASSERT_EQ(ordering.value()[1], &el1);
    ASSERT_EQ(el1, 1);
    ASSERT_EQ(el2, 2);
    ASSERT_TRUE(el2MovedLeft);
    ASSERT_TRUE(el1MovedRight);
}