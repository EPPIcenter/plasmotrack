//
// Created by Maxwell Murphy on 3/5/20.
//

#include "gtest/gtest.h"

#include "core/computation/OrderDerivedParentSet.h"
#include "core/parameters/Ordering.h"

using namespace transmission_nets::core::parameters;
using namespace transmission_nets::core::computation;

TEST(OrderingTest, HandlesSwapsNotifies) {
    auto el1 = std::make_shared<int>(1);
    auto el2 = std::make_shared<int>(2);
    auto el3 = std::make_shared<int>(3);
    auto el4 = std::make_shared<int>(4);

    std::unordered_map<std::shared_ptr<int>, int> movedLeftCounter{
            std::make_pair(el1, 0),
            std::make_pair(el2, 0),
            std::make_pair(el3, 0),
            std::make_pair(el4, 0),
    };

    std::unordered_map<std::shared_ptr<int>, int> movedRightCounter{
            std::make_pair(el1, 0),
            std::make_pair(el2, 0),
            std::make_pair(el3, 0),
            std::make_pair(el4, 0),
    };

    auto ordering = std::make_shared<Ordering<int>>();
    ordering->addElements({el1, el2, el3, el4});

    auto p_el1 = OrderDerivedParentSet<int, Ordering<int>>(ordering, el1, {el4, el2, el3});
    auto p_el2 = OrderDerivedParentSet<int, Ordering<int>>(ordering, el2, {el1, el3, el4});
    auto p_el3 = OrderDerivedParentSet<int, Ordering<int>>(ordering, el3, {el1, el2, el4});
    auto p_el4 = OrderDerivedParentSet<int, Ordering<int>>(ordering, el4, {el1, el2, el3});

    for (auto& el : ordering->value()) {
        ordering->add_keyed_moved_left_listener(el, [&]([[maybe_unused]] std::shared_ptr<int> tar) {
          movedLeftCounter[el] += 1;
        });
        ordering->add_keyed_moved_right_listener(el, [&]([[maybe_unused]] std::shared_ptr<int> tar) {
          movedRightCounter[el] += 1;
        });
    };


    ASSERT_EQ(ordering->value()[0], el1);
    ASSERT_EQ(ordering->value()[1], el2);
    ASSERT_EQ(*el1, 1);
    ASSERT_EQ(*el2, 2);

    ordering->saveState("state1");

    ordering->swap(0, 2);
    p_el1.serialize();
    p_el2.serialize();
    p_el3.serialize();
    p_el4.serialize();

    ordering->restoreState("state1");
    p_el1.serialize();
    p_el2.serialize();
    p_el3.serialize();
    p_el4.serialize();
}