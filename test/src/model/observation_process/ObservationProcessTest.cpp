//
// Created by Maxwell Murphy on 4/19/20.
//
//
#include "gtest/gtest.h"

#include "core/computation/Accumulator.h"

#include "core/containers/Infection.h"
#include "core/containers/Locus.h"

#include "core/datatypes/Alleles.h"

#include "core/parameters/Parameter.h"

#include "model/observation_process/AlleleCounter.h"
#include "model/observation_process/AlleleCounts.h"
#include "model/observation_process/ObservationProcessLikelihoodv1.h"

using namespace transmission_nets::core::computation;
using namespace transmission_nets::core::containers;
using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::parameters;
using namespace transmission_nets::model::observation_process;

constexpr int MAX_ALLELES = 32;

TEST(ObservationProcessTest, CoreTest) {
    using GeneticsImpl = AllelesBitSet<MAX_ALLELES>;
    using Infection = Infection<GeneticsImpl>;
    using AlleleCounter = AlleleCounter<GeneticsImpl>;
    using AlleleCounterAccumulator = Accumulator<AlleleCounter, AlleleCounts>;

    std::vector<std::shared_ptr<Locus>> loci{
            std::make_shared<Locus>("L1", 6),
            std::make_shared<Locus>("L2", 6),
            std::make_shared<Locus>("L3", 6)
    };

    std::vector<std::shared_ptr<Infection>> infections{};

    infections.reserve(4);
    for (int i = 0; i < 4; ++i) {
        auto infection = std::make_shared<Infection>(std::to_string(i), 1);
        infections.push_back(infection);
        for(auto &locus : loci) {
            infection->addGenetics(locus, "101010", "111111");
        }
    }

    std::vector<std::shared_ptr<AlleleCounter>> alleleCounters{};
    auto alleleCountAccumulator = std::make_shared<AlleleCounterAccumulator>();

    for (auto& infection : infections) {
        for (auto &[locus, obsGenotype] : infection->observedGenotype()) {
            alleleCounters.push_back(std::make_shared<AlleleCounter>(infection->latentGenotype(locus), obsGenotype));
            alleleCountAccumulator->addTarget(alleleCounters.back());
        }
    }

    EXPECT_EQ(alleleCountAccumulator->value().true_positive_count, 36);
    EXPECT_EQ(alleleCountAccumulator->value().true_negative_count, 0);
    EXPECT_EQ(alleleCountAccumulator->value().false_positive_count, 0);
    EXPECT_EQ(alleleCountAccumulator->value().false_negative_count, 36);


    auto falsePositiveRate = std::make_shared<Parameter<double>>(.01);
    auto falseNegativeRate = std::make_shared<Parameter<double>>(.05);

    ObservationProcessLikelihoodv1 target(alleleCountAccumulator, falsePositiveRate, falseNegativeRate);

    auto oldValue = target.value();

    infections.at(0)->latentGenotype(loci.at(0))->saveState("state1");
    infections.at(0)->latentGenotype(loci.at(0))->setValue(GeneticsImpl("100000"));

    EXPECT_EQ(alleleCountAccumulator->value().true_positive_count, 34);
    EXPECT_EQ(alleleCountAccumulator->value().true_negative_count, 3);
    EXPECT_EQ(alleleCountAccumulator->value().false_positive_count, 2);
    EXPECT_EQ(alleleCountAccumulator->value().false_negative_count, 33);

    EXPECT_TRUE(target.isDirty());
    auto newValue = target.value();
    EXPECT_FALSE(target.isDirty());

    EXPECT_GT(oldValue, newValue);

    infections.at(0)->latentGenotype(loci.at(0))->restoreState("state1");

    EXPECT_FALSE(target.isDirty());
    EXPECT_EQ(alleleCountAccumulator->value().true_positive_count, 36);
    EXPECT_EQ(alleleCountAccumulator->value().true_negative_count, 0);
    EXPECT_EQ(alleleCountAccumulator->value().false_positive_count, 0);
    EXPECT_EQ(alleleCountAccumulator->value().false_negative_count, 36);


}