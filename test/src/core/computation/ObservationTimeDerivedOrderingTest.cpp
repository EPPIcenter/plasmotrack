//
// Created by Maxwell Murphy on 4/7/21.
//


#include "core/computation/ObservationTimeDerivedOrdering.h"
#include "core/computation/OrderDerivedParentSet.h"
#include "core/containers/Infection.h"
#include "core/datatypes/Alleles.h"

#include "gtest/gtest.h"

using namespace transmission_nets::core::computation;
using namespace transmission_nets::core::containers;
using namespace transmission_nets::core::datatypes;

TEST(ObservationTimeDerivedOrderingTest, BaseTest) {
    using GeneticsImpl = AllelesBitSet<32>;
    using InfectionEventImpl = Infection<GeneticsImpl>;


    auto *inf1 = new InfectionEventImpl("1", 50);
    auto *inf2 = new InfectionEventImpl("2", 51);
    auto *inf3 = new InfectionEventImpl("3", 52);
    auto *inf4 = new InfectionEventImpl("4", 53);

    std::cout << "Making Ordering" << std::endl;
    auto *ord = new ObservationTimeDerivedOrdering<InfectionEventImpl>({inf1, inf2, inf3, inf4});

    OrderDerivedParentSet ps1(ord, inf1);
    OrderDerivedParentSet ps2(ord, inf2);
    OrderDerivedParentSet ps3(ord, inf3);
    OrderDerivedParentSet ps4(ord, inf4);

    ASSERT_EQ(ps1.value().size(), 0);
    ASSERT_EQ(ps2.value().size(), 1);
    ASSERT_EQ(ps3.value().size(), 2);
    ASSERT_EQ(ps4.value().size(), 3);

    inf1->infectionDuration().saveState("state1");
    inf1->infectionDuration().setValue(11);
    ASSERT_EQ(ps1.value().size(), 0);
    inf1->infectionDuration().restoreState("state1");
    ASSERT_EQ(ps1.value().size(), 0);

    inf2->infectionDuration().saveState("state1");
    inf2->infectionDuration().setValue(20);
    ASSERT_EQ(ps2.value().size(), 0);
    inf2->infectionDuration().restoreState("state1");
    ASSERT_EQ(ps2.value().size(), 1);

    inf3->infectionDuration().saveState("state1");
    inf3->infectionDuration().setValue(20);
    ASSERT_EQ(ps3.value().size(), 0);
    inf3->infectionDuration().restoreState("state1");
    ASSERT_EQ(ps3.value().size(), 2);

    inf4->infectionDuration().saveState("state1");
    inf4->infectionDuration().setValue(20);
    ASSERT_EQ(ps4.value().size(), 0);
    inf4->infectionDuration().restoreState("state1");
    ASSERT_EQ(ps4.value().size(), 3);

}