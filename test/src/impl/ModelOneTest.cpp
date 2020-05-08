//
// Created by Maxwell Murphy on 4/20/20.
//

#include "gtest/gtest.h"


#include "core/samplers/ConstrainedContinuousRandomWalk.h"
#include "core/samplers/SALTSampler.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler.h"
#include "core/samplers/OrderSampler.h"
#include "core/samplers/Scheduler.h"

#include "impl/model/ModelOne.h"
#include "impl/state/ModelOneState.h"


TEST(ModelOneTest, CoreTest) {
    using InfectionEvent = ModelOneState::InfectionEvent;
    using GeneticsImpl = ModelOneState::GeneticsImpl;

    using ZeroOneSampler = ConstrainedContinuousRandomWalk<0, 1, ModelOne, boost::random::mt19937>;
    using ZeroBoundedSampler = ConstrainedContinuousRandomWalk<0, std::numeric_limits<int>::max(), ModelOne, boost::random::mt19937>;
    using GeneticsSampler = RandomAllelesBitSetSampler<ModelOne, boost::random::mt19937, ModelOneState::GeneticsImpl>;

    ModelOneState state;

    for (int l = 1; l <= 12; ++l) {
        state.loci.push_back(new Locus("L" + std::to_string(l), 8));
    }

    for (int j = 0; j < 40; ++j) {
        state.infections.push_back(new InfectionEvent());
        for(auto &locus : state.loci) {
            state.infections.back()->addGenetics(locus, "10101010", "11111100");
        }
    }

    for(const auto &locus : state.loci) {
        state.alleleFrequencies.addLocus(locus);
    }

    state.infectionEventOrdering.addElements(state.infections);

    state.observationFalsePositiveRate.initializeValue(.025);
    state.observationFalseNegativeRate.initializeValue(.025);
    state.geometricGenerationProb.initializeValue(.5);
    state.ztMultiplicativeBinomialProb.initializeValue(.5);
    state.ztMultiplicativeBinomialAssoc.initializeValue(1.0);
    state.geometricCOIProb.initializeValue(.9);

    ModelOne model(state);

    state.observationFalsePositiveRate.saveState();
    auto currLik = model.value();
    state.observationFalsePositiveRate.setValue(.2);
    ASSERT_TRUE(model.isDirty());
    auto newLik = model.value();
    state.observationFalsePositiveRate.restoreState();
    auto oldLik = model.value();
    std::cout << "FPR Old/New: " << oldLik << " ||| " << newLik << std::endl;
    ASSERT_DOUBLE_EQ(currLik, oldLik);


    state.observationFalseNegativeRate.saveState();
    currLik = model.value();
    state.observationFalseNegativeRate.setValue(.2);
    ASSERT_TRUE(model.isDirty());
    newLik = model.value();
    ASSERT_FALSE(model.isDirty());
    state.observationFalseNegativeRate.restoreState();
    oldLik = model.value();
    std::cout << "FNR Old/New: " << oldLik << " ||| " << newLik << std::endl;
    ASSERT_DOUBLE_EQ(currLik, oldLik);

    state.geometricGenerationProb.saveState();
    currLik = model.value();
    state.geometricGenerationProb.setValue(.99999);
    ASSERT_TRUE(model.isDirty());
    newLik = model.value();
    ASSERT_FALSE(model.isDirty());
    state.geometricGenerationProb.restoreState();
    oldLik = model.value();
    std::cout << "Geo Gen Prob Old/New: " << oldLik << " ||| " << newLik << std::endl;
    ASSERT_DOUBLE_EQ(currLik, oldLik);

    state.ztMultiplicativeBinomialProb.saveState();
    currLik = model.value();
    state.ztMultiplicativeBinomialProb.setValue(.9);
    ASSERT_TRUE(model.isDirty());
    newLik = model.value();
    ASSERT_FALSE(model.isDirty());
    state.ztMultiplicativeBinomialProb.restoreState();
    oldLik = model.value();
    std::cout << "ZTMB Prob Old/New: " << oldLik << " ||| " << newLik << std::endl;
    ASSERT_DOUBLE_EQ(currLik, oldLik);

    state.ztMultiplicativeBinomialAssoc.saveState();
    currLik = model.value();
    state.ztMultiplicativeBinomialAssoc.setValue(.2);
    ASSERT_TRUE(model.isDirty());
    newLik = model.value();
    ASSERT_FALSE(model.isDirty());
    state.ztMultiplicativeBinomialAssoc.restoreState();
    oldLik = model.value();
    std::cout << "ZTMB Assoc Old/New: " << oldLik << " ||| " << newLik << std::endl;
    ASSERT_DOUBLE_EQ(currLik, oldLik);

    state.geometricCOIProb.saveState();
    currLik = model.value();
    state.geometricCOIProb.setValue(.99999);
    ASSERT_TRUE(model.isDirty());
    newLik = model.value();
    ASSERT_FALSE(model.isDirty());
    state.geometricCOIProb.restoreState();
    oldLik = model.value();
    std::cout << "Geo COI Prob Old/New: " << oldLik << " ||| " << newLik << std::endl;
    ASSERT_DOUBLE_EQ(currLik, oldLik);

    state.infections.at(10)->latentGenotype(state.loci.at(0)).saveState();
    currLik = model.value();
    state.infections.at(10)->latentGenotype(state.loci.at(0)).setValue(GeneticsImpl("00011100"));
    newLik = model.value();
    state.infections.at(10)->latentGenotype(state.loci.at(0)).restoreState();
    oldLik = model.value();
    std::cout << "Infection Old/New: " << oldLik << " ||| " << newLik << std::endl;
    ASSERT_DOUBLE_EQ(currLik, oldLik);


    boost::random::mt19937 r;
    Scheduler scheduler;
    scheduler.registerSampler(new OrderSampler(state.infectionEventOrdering, model, &r));
    scheduler.registerSampler(new ZeroOneSampler(state.observationFalsePositiveRate, model, &r));
    scheduler.registerSampler(new ZeroOneSampler(state.observationFalseNegativeRate, model, &r));
    scheduler.registerSampler(new ZeroOneSampler(state.geometricGenerationProb, model, &r));
    scheduler.registerSampler(new ZeroOneSampler(state.ztMultiplicativeBinomialProb, model, &r));
    scheduler.registerSampler(new ZeroBoundedSampler(state.ztMultiplicativeBinomialAssoc, model, &r));
    scheduler.registerSampler(new ZeroOneSampler(state.geometricCOIProb, model, &r));


    for(auto &infection : state.infections) {
        for(auto &locus : state.loci) {
            auto &latentGenotype = infection->latentGenotype(locus);
            scheduler.registerSampler(new GeneticsSampler(latentGenotype, model, &r));
        }
    }

    for(auto &locus : state.loci) {
        scheduler.registerSampler(new SALTSampler<ModelOne>(state.alleleFrequencies.alleleFrequencies(locus), model, &r));
    }

    for (int k = 0; k < 500; ++k) {
        scheduler.updateAndAdapt();
        std::cout << model.value() << std::endl;
    }
    std::cout << "Completed" << std::endl;
}