//
// Created by Maxwell Murphy on 3/6/20.
//

#include "core/parameters/Ordering.h"
#include "core/computation/OrderDerivedParentSet.h"
#include "gtest/gtest.h"

using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::computation;

TEST(OrderDerivedParentSetTest, HandlesReorder) {
    int el1 = 1;
    int el2 = 2;
    int el3 = 3;
    int el4 = 4;

    Ordering<int> ordering({&el1, &el2, &el3, &el4});
    OrderDerivedParentSet ops(ordering, el1);
    ASSERT_EQ(ops.value().size(), 0);

    ordering.saveState("state1");
    ordering.swap(0, 1);
    ASSERT_EQ(ops.value().size(), 1);
    ordering.acceptState();

    ordering.saveState("state1");
    ordering.swap(1, 2);
    ASSERT_EQ(ops.value().size(), 2);
    ordering.restoreState("state1");
    ASSERT_EQ(ops.value().size(), 1);

}