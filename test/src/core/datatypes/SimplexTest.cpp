//
// Created by Maxwell Murphy on 3/9/20.
//

#include <vector>

#include "gtest/gtest.h"
#include "core/datatypes/Simplex.h"
#include "core/parameters/Parameter.h"

using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::parameters;

TEST(SimplexTest, BasicTest) {
    Simplex av(3);
    ASSERT_DOUBLE_EQ(av.frequencies(0), 1.0/3.0);

    av.set({.33, .33, .33});
    ASSERT_DOUBLE_EQ(av.frequencies(0), 1.0 / 3.0);

    av.set({.25, .63, .33});

    Simplex av2({1, 2, 3});
    ASSERT_DOUBLE_EQ(av2.frequencies(0), 1.0 / 6.0);
}

TEST(SimplexTest, ParameterTest) {
    Simplex av(3);
    Parameter<Simplex> p(av);

    Parameter<Simplex> p2({.33, .33, .33});
    ASSERT_DOUBLE_EQ(p2.value().frequencies(0), 1.0/3.0);

    p2.saveState("state1");
    p2.setValue({.5, .6, .7});
    ASSERT_DOUBLE_EQ(p2.value().frequencies(0), .5 / (.5 + .6 + .7));
    p2.restoreState("state1");
    ASSERT_DOUBLE_EQ(p2.value().frequencies(0), 1.0/3.0);

    p2.saveState("state1");
    p2.setValue(Simplex({1.0, 2.0, 3.0}));
    ASSERT_DOUBLE_EQ(p2.value().frequencies(0), 1.0 / 6.0);
    p2.restoreState("state1");
    ASSERT_DOUBLE_EQ(p2.value().frequencies(0), 1.0/3.0);

    p2.saveState("state1");
    p2.setValue(Simplex(std::vector<double>({.1, .4, .6})));
    ASSERT_DOUBLE_EQ(p2.value().frequencies(0), .1 / (.1 + .4 + .6));
    p2.restoreState("state1");
    ASSERT_DOUBLE_EQ(p2.value().frequencies(0), 1.0/3.0);
}

