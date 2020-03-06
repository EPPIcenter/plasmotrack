//
// Created by Maxwell Murphy on 3/6/20.
//

#include "gtest/gtest.h"

#include "core/datatypes/Alleles.h"
#include "core/containers/Infection.h"


TEST(InfectionTest, HandlesChangedLatentAlleles) {
    using GeneticsImpl = AllelesBitSet<16>;

    Locus as1("AS1");
    bool allelesChanged = false;

    std::vector<LocusAssignment<GeneticsImpl, Locus>> dlas{{&as1, GeneticsImpl("011010")}};
    std::vector<LocusAssignment<GeneticsImpl, Locus>> plas{{&as1, GeneticsImpl("011010")}};
    Infection<GeneticsImpl> inf1(dlas, plas);

    inf1.add_post_change_listener([&]() { allelesChanged = true; });

    auto tmp = inf1.latentGenotype(&as1).value();
    tmp.flip(0);
    ASSERT_FALSE(allelesChanged);

    inf1.latentGenotype(&as1).saveState();
    inf1.latentGenotype(&as1).setValue(tmp);

    ASSERT_TRUE(allelesChanged);
    ASSERT_EQ(inf1.latentGenotype(&as1).value().allelesStr(), "011011");
    ASSERT_EQ(inf1.observedGenotype(&as1).value().allelesStr(), "011010");
}