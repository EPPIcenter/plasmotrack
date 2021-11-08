//
// Created by Maxwell Murphy on 3/6/20.
//

#include "gtest/gtest.h"

#include "core/datatypes/Alleles.h"
#include "core/containers/Infection.h"

using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::containers;
using namespace transmission_nets::core::parameters;

TEST(InfectionTest, HandlesChangedLatentAlleles) {
    using GeneticsImpl = AllelesBitSet<16>;
    using Infection = Infection<GeneticsImpl, Locus>;

    auto as1 = std::make_shared<Locus>("AS1", 6);
    auto as2 = std::make_shared<Locus>("AS2", 8);
    bool allelesChanged = false;

    auto inf1 = std::make_shared<Infection>("inf1", 10.0);
    inf1->addLatentGenetics(as1, GeneticsImpl("011010"));
    inf1->addLatentGenetics(as2, GeneticsImpl("00000011"));
    inf1->addObservedGenetics(as1, GeneticsImpl("011010"));
    inf1->addObservedGenetics(as2, GeneticsImpl("00000011"));
    inf1->add_post_change_listener([&]() {allelesChanged = true; });

    allelesChanged = false;
    auto tmp1 = inf1->latentGenotype(as1)->value();
    tmp1.flip(0);
    EXPECT_FALSE(allelesChanged);

    inf1->latentGenotype(as1)->saveState("State1");
    inf1->latentGenotype(as1)->setValue(tmp1);

    EXPECT_TRUE(allelesChanged);
    EXPECT_EQ(inf1->latentGenotype(as1)->value().allelesStr(), "111010");
    EXPECT_EQ(inf1->observedGenotype(as1)->value().allelesStr(), "011010");
    inf1->latentGenotype(as1)->acceptState();

    allelesChanged = false;
    auto tmp2 = inf1->latentGenotype(as2)->value();
    tmp2.flip(0);
    EXPECT_FALSE(allelesChanged);

    inf1->latentGenotype(as2)->saveState("State1");
    inf1->latentGenotype(as2)->setValue(tmp2);

    EXPECT_TRUE(allelesChanged);
    EXPECT_EQ(inf1->latentGenotype(as2)->value().allelesStr(), "10000011");
    EXPECT_EQ(inf1->observedGenotype(as2)->value().allelesStr(), "00000011");
    inf1->latentGenotype(as2)->acceptState();

    allelesChanged = false;
    inf1->latentGenotype(as2)->saveState("State1");
    inf1->latentGenotype(as2)->setValue(GeneticsImpl("00000011"));
    EXPECT_TRUE(allelesChanged);
    inf1->latentGenotype(as2)->acceptState();
}