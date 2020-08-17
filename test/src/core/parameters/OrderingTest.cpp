//
// Created by Maxwell Murphy on 3/5/20.
//

#include "gtest/gtest.h"

#include "core/parameters/Ordering.h"
#include "model/transmission_process/OrderDerivedParentSet.h"


TEST(OrderingTest, HandlesSwapsNotifies) {
    int el1 = 1;
    int el2 = 2;
    int el3 = 3;
    int el4 = 4;

    std::map<int*, int> movedLeftCounter{
            std::make_pair(&el1, 0),
            std::make_pair(&el2, 0),
            std::make_pair(&el3, 0),
            std::make_pair(&el4, 0),
    };

    std::map<int*, int> movedRightCounter{
            std::make_pair(&el1, 0),
            std::make_pair(&el2, 0),
            std::make_pair(&el3, 0),
            std::make_pair(&el4, 0),
    };

    Ordering<int> ordering;
    ordering.addElements({&el1, &el2, &el3, &el4});

    auto p_el1 = OrderDerivedParentSet<int>(ordering, el1);
    auto p_el2 = OrderDerivedParentSet<int>(ordering, el2);
    auto p_el3 = OrderDerivedParentSet<int>(ordering, el3);
    auto p_el4 = OrderDerivedParentSet<int>(ordering, el4);

    for(auto& el : ordering.value()) {
        ordering.add_keyed_moved_left_listener(el, [&]([[maybe_unused]] const int* tar) {
          movedLeftCounter[el] += 1;
        });
        ordering.add_keyed_moved_right_listener(el, [&]([[maybe_unused]] const int* tar) {
          movedRightCounter[el] += 1;
        });
    };


    ASSERT_EQ(ordering.value()[0], &el1);
    ASSERT_EQ(ordering.value()[1], &el2);
    ASSERT_EQ(el1, 1);
    ASSERT_EQ(el2, 2);

    ordering.saveState();

    ordering.swap(0, 2);
    p_el1.printSet();
    p_el2.printSet();
    p_el3.printSet();
    p_el4.printSet();

    ordering.restoreState();
    p_el1.printSet();
    p_el2.printSet();
    p_el3.printSet();
    p_el4.printSet();


}