//
// Created by Maxwell Murphy on 4/19/20.
//

#include "ModelOne.h"

ModelOne::ModelOne(ModelOneState state) : state(std::move(state)) {
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
            alleleCounters.push_back(std::make_unique<AlleleCounter>(infection->latentGenotype(locus), obsGenotype));
            alleleCountAccumulator.addTarget(*(alleleCounters.back()));
        }
        sourceTransmissionProcessList.push_back(std::make_unique<SourceTransmissionImpl>(
                *coiProb,
                state.alleleFrequencies,
                *infection
                                                )
        );
        parentSetList.push_back(new ParentSetImpl(state.infectionEventOrdering, *infection));
        transmissionProcessList.push_back(new TransmissionProcess(
                *nodeTransmissionProcess,
                *sourceTransmissionProcessList.back(),
                *infection,
                *parentSetList.back()
                                          )
        );
        likelihood.addTarget(*transmissionProcessList.back());
    }

    observationProcessLikelihood = new ObservationProcessLikelihood(
            alleleCountAccumulator,
            state.observationFalsePositiveRate,
            state.observationFalseNegativeRate
    );
    likelihood.addTarget(*observationProcessLikelihood);
}
