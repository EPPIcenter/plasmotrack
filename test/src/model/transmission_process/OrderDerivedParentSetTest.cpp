//
// Created by Maxwell Murphy on 3/6/20.
//

#include "core/computation/OrderDerivedParentSet.h"
#include "core/parameters/Ordering.h"
#include "gtest/gtest.h"

using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::computation;

TEST(OrderDerivedParentSetTest, HandlesReorder) {
//    auto el1 = std::make_shared<int>(1);
//    auto el2 = std::make_shared<int>(2);
//    auto el3 = std::make_shared<int>(3);
//    auto el4 = std::make_shared<int>(4);

    auto el1 = std::make_shared<Parameter<int>>(1);
    auto el2 = std::make_shared<Parameter<int>>(2);
    auto el3 = std::make_shared<Parameter<int>>(3);
    auto el4 = std::make_shared<Parameter<int>>(4);


    auto ordering = std::make_shared<Ordering<Parameter<int>>>(std::vector{el1, el2, el3, el4});
    auto ops      = std::make_shared<OrderDerivedParentSet<Parameter<int>, Ordering<Parameter<int>>>>(ordering, el1);
    ops->addAllowedParents({el2, el3, el4});
    ASSERT_EQ(ops->value().size(), 0);

    ordering->saveState(1);
    ordering->swap(0, 1);
    auto tmp = ops->value();
    ASSERT_EQ(ops->value().size(), 1);
    ordering->acceptState();

    ordering->saveState(1);
    ordering->swap(1, 2);
    ASSERT_EQ(ops->value().size(), 2);
    ordering->restoreState(1);
    ASSERT_EQ(ops->value().size(), 1);
}