//
// Created by Maxwell Murphy on 4/20/20.
//

#include "gtest/gtest.h"

#include "impl/model/ModelOne.h"
#include "impl/state/ModelOneState.h"


TEST(ModelOneTest, CoreTest) {
    using AlleleFrequencyContainer = ModelOneState::AlleleFrequencyContainer;
    using InfectionEvent = ModelOneState::InfectionEvent;
//    using GeneticsImpl = ModelOneState::GeneticsImpl;

    std::vector<Locus *> loci{
        new Locus("L1", 6),
        new Locus("L2", 6),
        new Locus("L3", 6),
        new Locus("L4", 6)
    };


    std::vector<InfectionEvent *> infections{};
    for (int j = 0; j < 10; ++j) {
        auto infection = new InfectionEvent();
        infections.push_back(infection);
        for(auto &locus : loci) {
            infection->addGenetics(locus, "101010", "111111");
        }
    }

    AlleleFrequencyContainer alleleFrequencies;
    for(const auto &locus : loci) {
        alleleFrequencies.addLocus(*locus);
    }

    Ordering<InfectionEvent> infectionEventOrdering;
    infectionEventOrdering.addElements(infections);

    Parameter<double> observationFalsePositiveRate(.05);
    Parameter<double> observationFalseNegativeRate(.05);

    Parameter<double> geometricGenerationProb(.5);
    Parameter<double> ztMultiplicativeBinomialProb(.5);
    Parameter<double> ztMultiplicativeBinomialAssoc(1.0);

    Parameter<double> geometricCOIProb(.9);


    ModelOneState state;

    state.loci = loci;
    state.infections = infections;
    state.alleleFrequencies = alleleFrequencies;
    state.infectionEventOrdering = infectionEventOrdering;
    state.observationFalsePositiveRate = observationFalsePositiveRate;
    state.observationFalseNegativeRate = observationFalseNegativeRate;
    state.geometricGenerationProb = geometricGenerationProb;
    state.ztMultiplicativeBinomialProb = ztMultiplicativeBinomialProb;
    state.ztMultiplicativeBinomialAssoc = ztMultiplicativeBinomialAssoc;
    state.geometricCOIProb = geometricCOIProb;

    ModelOne model(state);
    std::cout << model.value() << std::endl;


}