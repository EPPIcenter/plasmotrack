//
// Created by Maxwell Murphy on 4/19/20.
//

#include "ModelOne.h"


ModelOne::ModelOne(ModelOneState& state) : state(state) {
    coitp = new COITransitionProbImpl(state.ztMultiplicativeBinomialProb, state.ztMultiplicativeBinomialAssoc);
    intp = new InterTransmissionProbImpl(state.geometricGenerationProb);
    nodeTransmissionProcess = new NodeTransmissionImpl(*coitp, *intp);
    coiProb = new COIProbabilityImpl(state.geometricCOIProb);


    // Register Priors
    likelihood.addTarget(new BetaPrior(state.observationFalsePositiveRate, 2, 18));
    likelihood.addTarget(new BetaPrior(state.observationFalseNegativeRate, 2, 18));
    likelihood.addTarget(new GammaPrior(state.ztMultiplicativeBinomialAssoc, 10, .1));
    likelihood.addTarget(new BetaPrior(state.ztMultiplicativeBinomialProb, 80, 20));
    likelihood.addTarget(new BetaPrior(state.geometricCOIProb, 1, 1));
    likelihood.addTarget(new BetaPrior(state.geometricGenerationProb, 1, 1));

    for (auto &infection : state.infections) {
        for (auto &[locus, obsGenotype] : infection->observedGenotype()) {
            alleleCounters.push_back(new AlleleCounter(infection->latentGenotype(locus), obsGenotype));
            alleleCountAccumulator.addTarget(alleleCounters.back());
        }
        parentSetList.push_back(new ParentSetImpl(state.infectionEventOrdering, *infection));
    }

    observationProcessLikelihood = new ObservationProcessLikelihood(
            alleleCountAccumulator,
            state.observationFalsePositiveRate,
            state.observationFalseNegativeRate
    );
    likelihood.addTarget(observationProcessLikelihood);

    for (unsigned int j = 0; j < state.infections.size(); ++j) {
        auto infection = state.infections[j];
        auto parentSet = parentSetList[j];

        sourceTransmissionProcessList.push_back(new SourceTransmissionImpl(
                *coiProb,
                state.alleleFrequencies,
                *infection
                                                )
        );

        transmissionProcessList.push_back(new TransmissionProcess(
                *nodeTransmissionProcess,
                *sourceTransmissionProcessList.back(),
                *infection,
                *parentSet
                                          )
        );

        likelihood.addTarget(transmissionProcessList.back());

    }
}

double ModelOne::value() {
    return likelihood.value();
}

bool ModelOne::isDirty() {
    return likelihood.isDirty();
}
