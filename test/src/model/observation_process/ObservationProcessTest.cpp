//
// Created by Maxwell Murphy on 4/19/20.
//

#include "gtest/gtest.h"

#include "core/computation/Accumulator.h"

#include "core/containers/Infection.h"
#include "core/containers/Locus.h"

#include "core/datatypes/Alleles.h"

#include "model/observation_process/AlleleCounts.h"
#include "model/observation_process/AlleleCounter.h"
#include "model/observation_process/ObservationProcessLikelihood.h"


constexpr int MAX_ALLELES = 32;

TEST(ObservationProcessTest, CoreTest) {
    using GeneticsImpl = AllelesBitSet<MAX_ALLELES>;
    using Infection = Infection<GeneticsImpl>;
    using AlleleCounter = AlleleCounter<GeneticsImpl>;
    using AlleleCounterAccumulator = Accumulator<AlleleCounter, AlleleCounts>;

    std::vector<Locus *> loci{
            new Locus("L1", 6),
            new Locus("L2", 6),
            new Locus("L3", 6)
    };

    std::vector<Infection *> infections{};

    infections.reserve(4);
    for (int i = 0; i < 4; ++i) {
        auto infection = new Infection(std::to_string(i));
        infections.push_back(infection);
        for(auto &locus : loci) {
            infection->addGenetics(locus, "101010", "111111");
        }
    }

    std::vector<AlleleCounter *> alleleCounters{};
    AlleleCounterAccumulator alleleCountAccumulator;

    for (auto& infection : infections) {
        for (auto &[locus, obsGenotype] : infection->observedGenotype()) {
            alleleCounters.push_back(new AlleleCounter(infection->latentGenotype(locus), obsGenotype));
            alleleCountAccumulator.addTarget(alleleCounters.back());
        }
    }

    EXPECT_EQ(alleleCountAccumulator.value().true_positive_count, 36);
    EXPECT_EQ(alleleCountAccumulator.value().true_negative_count, 0);
    EXPECT_EQ(alleleCountAccumulator.value().false_positive_count, 0);
    EXPECT_EQ(alleleCountAccumulator.value().false_negative_count, 36);


    Parameter<double> falsePositiveRate(.01);
    Parameter<double> falseNegativeRate(.05);

    ObservationProcessLikelihood target(alleleCountAccumulator, falsePositiveRate, falseNegativeRate);

    auto oldValue = target.value();

    infections.at(0)->latentGenotype(loci.at(0)).saveState("state1");
    infections.at(0)->latentGenotype(loci.at(0)).setValue(GeneticsImpl("100000"));

    EXPECT_EQ(alleleCountAccumulator.value().true_positive_count, 34);
    EXPECT_EQ(alleleCountAccumulator.value().true_negative_count, 3);
    EXPECT_EQ(alleleCountAccumulator.value().false_positive_count, 2);
    EXPECT_EQ(alleleCountAccumulator.value().false_negative_count, 33);

    EXPECT_TRUE(target.isDirty());
    auto newValue = target.value();
    EXPECT_FALSE(target.isDirty());

    EXPECT_GT(oldValue, newValue);

    infections.at(0)->latentGenotype(loci.at(0)).restoreState("state1");

    EXPECT_FALSE(target.isDirty());
    EXPECT_EQ(alleleCountAccumulator.value().true_positive_count, 36);
    EXPECT_EQ(alleleCountAccumulator.value().true_negative_count, 0);
    EXPECT_EQ(alleleCountAccumulator.value().false_positive_count, 0);
    EXPECT_EQ(alleleCountAccumulator.value().false_negative_count, 36);


}