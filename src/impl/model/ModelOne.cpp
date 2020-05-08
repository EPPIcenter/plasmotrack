//
// Created by Maxwell Murphy on 4/19/20.
//

#include "ModelOne.h"


ModelOne::ModelOne(ModelOneState& state) : state(state) {
    init();
}

double ModelOne::value() {
    return likelihood.value();
}

void ModelOne::init() {
    coitp = new COITransitionProbImpl(state.ztMultiplicativeBinomialProb, state.ztMultiplicativeBinomialAssoc);
    intp = new InterTransmissionProbImpl(state.geometricGenerationProb);
    nodeTransmissionProcess = new NodeTransmissionImpl(*coitp, *intp);
    coiProb = new COIProbabilityImpl(state.geometricCOIProb);

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

        auto logTransmissionProcess = new LogTransformer(*(transmissionProcessList.back()));
        logTransmissionProcessList.push_back(logTransmissionProcess);
        likelihood.addTarget(logTransmissionProcessList.back());

        likelihood.addTarget(new BetaPrior(state.observationFalsePositiveRate, 5, 95));
        likelihood.addTarget(new BetaPrior(state.observationFalseNegativeRate, 5, 95));
    }

}

bool ModelOne::isDirty() {
    return likelihood.isDirty();
}
