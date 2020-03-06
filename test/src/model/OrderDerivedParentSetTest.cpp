//
// Created by Maxwell Murphy on 3/6/20.
//

#include "gtest/gtest.h"
#include "model/transmission_process/OrderDerivedParentSet.h"

TEST(OrderDerivedParentSetTest, HandlesReorder) {
    auto el1 = 1;
    auto el2 = 2;
    auto el3 = 3;
    auto el4 = 4;

    Ordering<int> ordering({&el1, &el2, &el3, &el4});
    OrderDerivedParentSet ops(ordering, el1);
    ASSERT_EQ(ops.value().size(), 0);

    ordering.saveState();
    ordering.swap(0, 1);
    ASSERT_EQ(ops.value().size(), 1);
    ordering.acceptState();

    ordering.saveState();
    ordering.swap(1, 2);
    ASSERT_EQ(ops.value().size(), 2);
    ordering.restoreState();
    ASSERT_EQ(ops.value().size(), 1);

}