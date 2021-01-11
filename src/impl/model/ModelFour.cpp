//
// Created by Maxwell Murphy on 12/10/20.
//

#include "ModelFour.h"

#include "core/distributions/pdfs/BetaLogPDF.h"
#include "core/distributions/pdfs/GammaLogPDF.h"


namespace transmission_nets::impl {

    ModelFour::ModelFour(std::map<std::string, LocusImpl *>& loci,
                       std::vector<InfectionEvent *>& infections,
                       std::map<InfectionEvent *, std::vector<InfectionEvent *>>& disallowedParents) : state(loci, infections, disallowedParents) {
        intp = new InterTransmissionProbImpl(state.geometricGenerationProb);
        nodeTransmissionProcess = new NodeTransmissionImpl(state.mutationProb, state.lossProb, *intp);
        coiProb = new COIProbabilityImpl(state.meanCOI);

        // Register Priors
//        likelihood.addTarget(new core::distributions::BetaLogPDF(state.observationFalsePositiveRate, 1, 100));
//        likelihood.addTarget(new core::distributions::BetaLogPDF(state.observationFalseNegativeRate, 1, 100));
        likelihood.addTarget(new core::distributions::BetaLogPDF(state.mutationProb, 1, 500));
        likelihood.addTarget(new core::distributions::BetaLogPDF(state.lossProb, 1, 1));
        likelihood.addTarget(new core::distributions::GammaLogPDF(state.meanCOI, 1, 5));
        likelihood.addTarget(new core::distributions::BetaLogPDF(state.geometricGenerationProb, 1, 1));
        for (auto& obs : state.observationFalsePositiveRates) {
            likelihood.addTarget(new core::distributions::BetaLogPDF(obs, 1, 99));
        }
        for (auto& obs : state.observationFalseNegativeRates) {
            likelihood.addTarget(new core::distributions::BetaLogPDF(obs, 10, 90));
        }


        int i = 0;
        for (auto &infection : state.infections) {
            alleleCountAccumulators.push_back(new AlleleCounterAccumulator());
            for (auto &[locus, obsGenotype] : infection->observedGenotype()) {
                alleleCounters.push_back(new AlleleCounterImpl(infection->latentGenotype(locus), obsGenotype));
                alleleCountAccumulators.back()->addTarget(alleleCounters.back());
            }
            observationProcessLikelihood = new model::observation_process::ObservationProcessLikelihood(
                    *(alleleCountAccumulators.back()),
                    state.observationFalseNegativeRates[i],
                    state.observationFalsePositiveRates[i]);
            likelihood.addTarget(observationProcessLikelihood);
            parentSetList.push_back(new ParentSetImpl(state.infectionEventOrdering, *infection));
            parentSetList.back()->addDisallowedParents(state.disallowedParents.at(infection));
            i++;
        }

        for (unsigned int j = 0; j < state.infections.size(); ++j) {
            auto infection = state.infections[j];
            auto parentSet = parentSetList[j];

            sourceTransmissionProcessList.push_back(new SourceTransmissionImpl(
                    *coiProb,
                    state.alleleFrequencies,
                    *infection));

            transmissionProcessList.push_back(new TransmissionProcess(
                    *nodeTransmissionProcess,
                    *sourceTransmissionProcessList.back(),
                    *infection,
                    *parentSet));

            likelihood.addTarget(transmissionProcessList.back());
        }
    }

    ModelFour::Likelihood ModelFour::value() {
        return likelihood.value();
    }

    bool ModelFour::isDirty() {
        return likelihood.isDirty();
    }


    ModelFour::State::State(std::map<std::string, LocusImpl *>& loci,
                           std::vector<InfectionEvent *>& infections,
                           std::map<InfectionEvent *, std::vector<InfectionEvent *>>& disallowedParents) : loci(loci), infections(infections), disallowedParents(disallowedParents) {
        for (const auto &[locus_label, locus] : loci) {
            alleleFrequencies.addLocus(locus);
        }

        infectionEventOrdering.addElements(infections);

        for (size_t _ = 0; _ < infections.size(); ++_) {
            observationFalsePositiveRates.emplace_back();
            observationFalsePositiveRates.back().initializeValue(.01);
            observationFalseNegativeRates.emplace_back();
            observationFalseNegativeRates.back().initializeValue(.1);
        }

//        observationFalsePositiveRate.initializeValue(.001);
//        observationFalseNegativeRate.initializeValue(.1);

        geometricGenerationProb.initializeValue(.9);
        lossProb.initializeValue(.1);
        mutationProb.initializeValue(1e-8);
        meanCOI.initializeValue(5);
    }
}