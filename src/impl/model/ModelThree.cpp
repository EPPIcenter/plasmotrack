//
// Created by Maxwell Murphy on 6/18/20.
//

#include "ModelThree.h"

ModelThree::ModelThree(ModelThreeState& state) : state(state) {
    intp = new InterTransmissionProbImpl(state.geometricGenerationProb);
    nodeTransmissionProcess = new NodeTransmissionImpl(state.mutationProb, state.lossProb, *intp);
    coiProb = new COIProbabilityImpl(state.meanCOI);


    // Register Priors
    likelihood.addTarget(new BetaPrior(state.observationFalsePositiveRate, 1, 1));
    likelihood.addTarget(new BetaPrior(state.observationFalseNegativeRate, 1, 1));
    likelihood.addTarget(new BetaPrior(state.mutationProb, 1, 1));
    likelihood.addTarget(new BetaPrior(state.lossProb, 1, 1));
    likelihood.addTarget(new GammaPrior(state.meanCOI, 1, 1));
    likelihood.addTarget(new BetaPrior(state.geometricGenerationProb, 1, 1));

    for (auto &infection : state.infections) {
        for (auto &[locus, obsGenotype] : infection->observedGenotype()) {
            alleleCounters.push_back(new AlleleCounter(infection->latentGenotype(locus), obsGenotype));
            alleleCountAccumulator.addTarget(alleleCounters.back());
        }
    }

    observationProcessLikelihood = new ObservationProcessLikelihood(
            alleleCountAccumulator,
            state.observationFalsePositiveRate,
            state.observationFalseNegativeRate
    );
    likelihood.addTarget(observationProcessLikelihood);

    for (unsigned int j = 0; j < state.infections.size(); ++j) {
        auto infection = state.infections[j];
        auto parentSet = state.transmissionNetwork.parentSet(infection);

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

double ModelThree::value() {
    return likelihood.value();
}

bool ModelThree::isDirty() {
    return likelihood.isDirty();
}