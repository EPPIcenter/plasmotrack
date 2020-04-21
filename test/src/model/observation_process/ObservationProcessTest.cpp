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

//    std::vector<std::shared_ptr<Locus>> loci{
//        std::make_shared<Locus>("AS1", 5),
//        std::make_shared<Locus>("AS2", 6)
//    };

    Locus as1("AS1", 5);
    Locus as2("AS2", 6);

    std::vector<Infection::LocusGeneticsAssignment> observed{
            {&as1, GeneticsImpl("11010")},
            {&as2, GeneticsImpl("000011")}
    };

    std::vector<Infection::LocusGeneticsAssignment> latent{
            {&as1, GeneticsImpl("11010")},
            {&as2, GeneticsImpl("000011")}
    };

    std::vector<std::unique_ptr<Infection>> infections{};

    infections.reserve(4);
    for (int i = 0; i < 4; ++i) {
        infections.push_back(std::make_unique<Infection>(latent, observed));
    }

    std::vector<std::unique_ptr<AlleleCounter>> alleleCounters{};
    AlleleCounterAccumulator alleleCountAccumulator;

    for (auto& infection : infections) {
        for (auto &[locus, obsGenotype] : infection->observedGenotype()) {
            alleleCounters.push_back(std::make_unique<AlleleCounter>(infection->latentGenotype(locus), obsGenotype));
            alleleCountAccumulator.addTarget(*(alleleCounters.back()));
        }
    }

    EXPECT_EQ(alleleCountAccumulator.value().true_positive_count, 20);
    EXPECT_EQ(alleleCountAccumulator.value().true_negative_count, 24);
    EXPECT_EQ(alleleCountAccumulator.value().false_positive_count, 0);
    EXPECT_EQ(alleleCountAccumulator.value().false_negative_count, 0);


    Parameter<double> falsePositiveRate(.01);
    Parameter<double> falseNegativeRate(.05);

    ObservationProcessLikelihood target(alleleCountAccumulator, falsePositiveRate, falseNegativeRate);

    auto oldValue = target.value();

    infections.at(0)->latentGenotype(&as1).saveState();
    infections.at(0)->latentGenotype(&as1).setValue(GeneticsImpl("11111"));

    EXPECT_EQ(alleleCountAccumulator.value().true_positive_count, 20);
    EXPECT_EQ(alleleCountAccumulator.value().true_negative_count, 22);
    EXPECT_EQ(alleleCountAccumulator.value().false_positive_count, 0);
    EXPECT_EQ(alleleCountAccumulator.value().false_negative_count, 2);

    EXPECT_TRUE(target.isDirty());
    auto newValue = target.value();
    EXPECT_FALSE(target.isDirty());

    EXPECT_GT(oldValue, newValue);

    infections.at(0)->latentGenotype(&as1).restoreState();

    EXPECT_FALSE(target.isDirty());
    EXPECT_EQ(alleleCountAccumulator.value().true_positive_count, 20);
    EXPECT_EQ(alleleCountAccumulator.value().true_negative_count, 24);
    EXPECT_EQ(alleleCountAccumulator.value().false_positive_count, 0);
    EXPECT_EQ(alleleCountAccumulator.value().false_negative_count, 0);


}