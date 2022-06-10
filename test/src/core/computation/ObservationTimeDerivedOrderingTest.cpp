//
// Created by Maxwell Murphy on 4/7/21.
//


#include "core/computation/ObservationTimeDerivedOrdering.h"
#include "core/computation/OrderDerivedParentSet.h"
#include "core/containers/Infection.h"
#include "core/datatypes/Alleles.h"
#include "core/io/serialize.h"

#include "gtest/gtest.h"

#include <memory>
#include <vector>

using namespace transmission_nets::core::computation;
using namespace transmission_nets::core::containers;
using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::io;

TEST(ObservationTimeDerivedOrderingTest, BaseTest) {
    using GeneticsImpl       = AllelesBitSet<32>;
    using InfectionEventImpl = Infection<GeneticsImpl>;


    auto inf1       = std::make_shared<InfectionEventImpl>("1", 1000);
    auto inf2       = std::make_shared<InfectionEventImpl>("2", 1010);
    auto inf3       = std::make_shared<InfectionEventImpl>("3", 1020);
    auto inf4       = std::make_shared<InfectionEventImpl>("4", 1030);
    auto infections = std::vector({inf1, inf2, inf3, inf4});
    auto ord        = std::make_shared<ObservationTimeDerivedOrdering<InfectionEventImpl>>(infections);

    auto inf1_parents = std::vector({inf2, inf3, inf4});
    auto inf2_parents = std::vector({inf1, inf3, inf4});
    auto inf3_parents = std::vector({inf1, inf2, inf4});
    auto inf4_parents = std::vector({inf1, inf2, inf3});

    bool ps1_parent_added = false;
    bool ps2_parent_added = false;
    bool ps3_parent_added = false;
    bool ps4_parent_added = false;

    bool ps1_parent_removed = false;
    bool ps2_parent_removed = false;
    bool ps3_parent_removed = false;
    bool ps4_parent_removed = false;

    OrderDerivedParentSet ps1(ord, inf1, inf1_parents);
    ps1.add_element_added_listener([&]([[maybe_unused]] auto el) { ps1_parent_added = true; });
    ps1.add_element_removed_listener([&]([[maybe_unused]] auto el) { ps1_parent_removed = true; });

    OrderDerivedParentSet ps2(ord, inf2, inf2_parents);
    ps2.add_element_added_listener([&]([[maybe_unused]] auto el) { ps2_parent_added = true; });
    ps2.add_element_removed_listener([&]([[maybe_unused]] auto el) { ps2_parent_removed = true; });

    OrderDerivedParentSet ps3(ord, inf3, inf3_parents);
    ps3.add_element_added_listener([&]([[maybe_unused]] auto el) { ps3_parent_added = true; });
    ps3.add_element_removed_listener([&]([[maybe_unused]] auto el) { ps3_parent_removed = true; });

    OrderDerivedParentSet ps4(ord, inf4, inf4_parents);
    ps4.add_element_added_listener([&]([[maybe_unused]] auto el) { ps4_parent_added = true; });
    ps4.add_element_removed_listener([&]([[maybe_unused]] auto el) { ps4_parent_removed = true; });

    ASSERT_EQ(ps1.value().size(), 0);
    ASSERT_EQ(ps2.value().size(), 1);
    ASSERT_EQ(ps3.value().size(), 2);
    ASSERT_EQ(ps4.value().size(), 3);

    inf1->infectionDuration()->saveState("state1");
    inf1->infectionDuration()->setValue(inf1->infectionDuration()->value() - 5);
    ASSERT_EQ(ps1.value().size(), 0);
    ASSERT_FALSE(ps1_parent_added);
    ASSERT_FALSE(ps2_parent_added);
    ASSERT_FALSE(ps3_parent_added);
    ASSERT_FALSE(ps4_parent_added);
    ASSERT_FALSE(ps1_parent_removed);
    ASSERT_FALSE(ps2_parent_removed);
    ASSERT_FALSE(ps3_parent_removed);
    ASSERT_FALSE(ps4_parent_removed);
    inf1->infectionDuration()->restoreState("state1");
    ASSERT_EQ(ps1.value().size(), 0);

    inf1->infectionDuration()->saveState("state1");
    inf1->infectionDuration()->setValue(inf1->infectionDuration()->value() - 15);
    ASSERT_EQ(ps1.value().size(), 1);
    ASSERT_TRUE(ps1_parent_added);
    ASSERT_FALSE(ps2_parent_added);
    ASSERT_FALSE(ps3_parent_added);
    ASSERT_FALSE(ps4_parent_added);
    ASSERT_FALSE(ps1_parent_removed);
    ASSERT_TRUE(ps2_parent_removed);
    ASSERT_FALSE(ps3_parent_removed);
    ASSERT_FALSE(ps4_parent_removed);
    inf1->infectionDuration()->restoreState("state1");
    ASSERT_EQ(ps1.value().size(), 0);
    ps1_parent_added   = false;
    ps2_parent_added   = false;
    ps3_parent_added   = false;
    ps4_parent_added   = false;
    ps1_parent_removed = false;
    ps2_parent_removed = false;
    ps3_parent_removed = false;
    ps4_parent_removed = false;


    inf1->infectionDuration()->saveState("state1");
    inf1->infectionDuration()->setValue(inf1->infectionDuration()->value() - 25);
    ASSERT_EQ(ps1.value().size(), 2);
    ASSERT_TRUE(ps1_parent_added);
    ASSERT_FALSE(ps2_parent_added);
    ASSERT_FALSE(ps3_parent_added);
    ASSERT_FALSE(ps4_parent_added);
    ASSERT_FALSE(ps1_parent_removed);
    ASSERT_TRUE(ps2_parent_removed);
    ASSERT_TRUE(ps3_parent_removed);
    ASSERT_FALSE(ps4_parent_removed);
    inf1->infectionDuration()->restoreState("state1");
    ASSERT_EQ(ps1.value().size(), 0);
    ps1_parent_added   = false;
    ps2_parent_added   = false;
    ps3_parent_added   = false;
    ps4_parent_added   = false;
    ps1_parent_removed = false;
    ps2_parent_removed = false;
    ps3_parent_removed = false;
    ps4_parent_removed = false;

    inf1->infectionDuration()->saveState("state1");
    inf1->infectionDuration()->setValue(inf1->infectionDuration()->value() - 35);
    ASSERT_EQ(ps1.value().size(), 3);
    inf1->infectionDuration()->restoreState("state1");
    ASSERT_EQ(ps1.value().size(), 0);

    inf2->infectionDuration()->saveState("state1");
    inf2->infectionDuration()->setValue(inf2->infectionDuration()->value() + 15);
    ASSERT_EQ(ps2.value().size(), 0);
    inf2->infectionDuration()->restoreState("state1");
    ASSERT_EQ(ps2.value().size(), 1);

    inf2->infectionDuration()->saveState("state1");
    inf2->infectionDuration()->setValue(inf2->infectionDuration()->value() + 5);
    ASSERT_EQ(ps2.value().size(), 1);
    inf2->infectionDuration()->restoreState("state1");
    ASSERT_EQ(ps2.value().size(), 1);

    inf2->infectionDuration()->saveState("state1");
    inf2->infectionDuration()->setValue(inf2->infectionDuration()->value() - 15);
    ASSERT_EQ(ps2.value().size(), 2);
    inf2->infectionDuration()->restoreState("state1");
    ASSERT_EQ(ps2.value().size(), 1);

    inf2->infectionDuration()->saveState("state1");
    inf2->infectionDuration()->setValue(inf2->infectionDuration()->value() - 25);
    ASSERT_EQ(ps2.value().size(), 3);
    inf2->infectionDuration()->restoreState("state1");
    ASSERT_EQ(ps2.value().size(), 1);


    inf3->infectionDuration()->saveState("state1");
    inf3->infectionDuration()->setValue(inf3->infectionDuration()->value() + 25);
    ASSERT_EQ(ps3.value().size(), 0);
    inf3->infectionDuration()->restoreState("state1");
    ASSERT_EQ(ps3.value().size(), 2);

    inf3->infectionDuration()->saveState("state1");
    inf3->infectionDuration()->setValue(inf3->infectionDuration()->value() + 15);
    ASSERT_EQ(ps3.value().size(), 1);
    inf3->infectionDuration()->restoreState("state1");
    ASSERT_EQ(ps3.value().size(), 2);

    inf3->infectionDuration()->saveState("state1");
    inf3->infectionDuration()->setValue(inf3->infectionDuration()->value() + 5);
    ASSERT_EQ(ps3.value().size(), 2);
    inf3->infectionDuration()->restoreState("state1");
    ASSERT_EQ(ps3.value().size(), 2);

    inf3->infectionDuration()->saveState("state1");
    inf3->infectionDuration()->setValue(inf3->infectionDuration()->value() - 5);
    ASSERT_EQ(ps3.value().size(), 2);
    inf3->infectionDuration()->restoreState("state1");
    ASSERT_EQ(ps3.value().size(), 2);

    inf3->infectionDuration()->saveState("state1");
    inf3->infectionDuration()->setValue(inf3->infectionDuration()->value() - 15);
    ASSERT_EQ(ps3.value().size(), 3);
    inf3->infectionDuration()->restoreState("state1");
    ASSERT_EQ(ps3.value().size(), 2);
    ps1_parent_added   = false;
    ps2_parent_added   = false;
    ps3_parent_added   = false;
    ps4_parent_added   = false;
    ps1_parent_removed = false;
    ps2_parent_removed = false;
    ps3_parent_removed = false;
    ps4_parent_removed = false;


    inf4->infectionDuration()->saveState("state1");
    inf4->infectionDuration()->setValue(inf4->infectionDuration()->value() + 35);
    ASSERT_EQ(ps4.value().size(), 0);
    ASSERT_TRUE(ps1_parent_added);
    ASSERT_TRUE(ps2_parent_added);
    ASSERT_TRUE(ps3_parent_added);
    ASSERT_FALSE(ps4_parent_added);
    ASSERT_FALSE(ps1_parent_removed);
    ASSERT_FALSE(ps2_parent_removed);
    ASSERT_FALSE(ps3_parent_removed);
    ASSERT_TRUE(ps4_parent_removed);
    inf4->infectionDuration()->restoreState("state1");
    ASSERT_EQ(ps4.value().size(), 3);

    inf4->infectionDuration()->saveState("state1");
    inf4->infectionDuration()->setValue(inf4->infectionDuration()->value() + 25);
    ASSERT_EQ(ps4.value().size(), 1);
    inf4->infectionDuration()->restoreState("state1");
    ASSERT_EQ(ps4.value().size(), 3);

    inf4->infectionDuration()->saveState("state1");
    inf4->infectionDuration()->setValue(inf4->infectionDuration()->value() + 15);
    ASSERT_EQ(ps4.value().size(), 2);
    inf4->infectionDuration()->restoreState("state1");
    ASSERT_EQ(ps4.value().size(), 3);

    inf4->infectionDuration()->saveState("state1");
    inf4->infectionDuration()->setValue(inf4->infectionDuration()->value() - 5);
    ASSERT_EQ(ps4.value().size(), 3);
    inf4->infectionDuration()->restoreState("state1");
    ASSERT_EQ(ps4.value().size(), 3);

    inf4->infectionDuration()->saveState("state1");
    inf4->infectionDuration()->setValue(inf4->infectionDuration()->value() - 15);
    ASSERT_EQ(ps4.value().size(), 3);
    inf4->infectionDuration()->restoreState("state1");
    ASSERT_EQ(ps4.value().size(), 3);
}