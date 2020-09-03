//
// Created by Maxwell Murphy on 6/5/20.
//

#include "ModelTwo.h"

#include "core/distributions/pdfs/BetaLogPDF.h"
#include "core/distributions/pdfs/GammaLogPDF.h"


ModelTwo::ModelTwo(ModelTwoState& state) : state(state) {
    intp = new InterTransmissionProbImpl(state.geometricGenerationProb);
    nodeTransmissionProcess = new NodeTransmissionImpl(state.mutationProb, state.lossProb, *intp);
    coiProb = new COIProbabilityImpl(state.meanCOI);


    // Register Priors
    likelihood.addTarget(new BetaLogPDF(state.observationFalsePositiveRate, 1, 100));
    likelihood.addTarget(new BetaLogPDF(state.observationFalseNegativeRate, 1, 100));
    likelihood.addTarget(new BetaLogPDF(state.mutationProb, 1, 500));
    likelihood.addTarget(new BetaLogPDF(state.lossProb, 1, 1));
    likelihood.addTarget(new GammaLogPDF(state.meanCOI, 1, 1));
    likelihood.addTarget(new BetaLogPDF(state.geometricGenerationProb, 1, 1));

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

double ModelTwo::value() {
    return likelihood.value();
}

bool ModelTwo::isDirty() {
    return likelihood.isDirty();
}
