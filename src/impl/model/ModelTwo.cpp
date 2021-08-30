////
//// Created by Maxwell Murphy on 6/5/20.
////
//
//#include "ModelTwo.h"
//
//#include <utility>
//
//#include "core/distributions/pdfs/BetaLogPDF.h"
//#include "core/distributions/pdfs/GammaLogPDF.h"
//
//
//namespace transmission_nets::impl {
//
//    ModelTwo::ModelTwo(std::map<std::string, LocusImpl *>& loci,
//                       std::vector<InfectionEvent *>& infections,
//                       std::map<InfectionEvent *, std::vector<InfectionEvent *>>& disallowedParents) : state(loci, infections, disallowedParents) {
//        intp = new InterTransmissionProbImpl(state.geometricGenerationProb);
//        nodeTransmissionProcess = new NodeTransmissionImpl(state.mutationProb, state.lossProb, *intp);
//        coiProb = new COIProbabilityImpl(state.meanCOI);
//
//        // Register Priors
//        likelihood.addTarget(new core::distributions::BetaLogPDF(state.observationFalsePositiveRate, 1, 100));
//        likelihood.addTarget(new core::distributions::BetaLogPDF(state.observationFalseNegativeRate, 1, 100));
//        likelihood.addTarget(new core::distributions::BetaLogPDF(state.mutationProb, 1, 500));
//        likelihood.addTarget(new core::distributions::BetaLogPDF(state.lossProb, 1, 1));
//        likelihood.addTarget(new core::distributions::GammaLogPDF(state.meanCOI, 1, 1));
//        likelihood.addTarget(new core::distributions::BetaLogPDF(state.geometricGenerationProb, 1, 1));
//
//        for (auto &infection : state.infections) {
//            for (auto &[locus, obsGenotype] : infection->observedGenotype()) {
//                alleleCounters.push_back(new AlleleCounterImpl(infection->latentGenotype(locus), obsGenotype));
//                alleleCountAccumulator.addTarget(alleleCounters.back());
//            }
//            parentSetList.push_back(new ParentSetImpl(&(state.infectionEventOrdering), infection));
//            parentSetList.back()->addDisallowedParents(state.disallowedParents.at(infection));
//        }
//
//        observationProcessLikelihood = new model::observation_process::ObservationProcessLikelihood(
//                alleleCountAccumulator,
//                state.observationFalsePositiveRate,
//                state.observationFalseNegativeRate);
//        likelihood.addTarget(observationProcessLikelihood);
//
//        for (unsigned int j = 0; j < state.infections.size(); ++j) {
//            auto infection = state.infections[j];
//            auto parentSet = parentSetList[j];
//
//            sourceTransmissionProcessList.push_back(new SourceTransmissionImpl(
//                    *coiProb,
//                    state.alleleFrequencies,
//                    *infection));
//
//            transmissionProcessList.push_back(new TransmissionProcess(
//                    *nodeTransmissionProcess,
//                    *sourceTransmissionProcessList.back(),
//                    *infection,
//                    *parentSet));
//
//            likelihood.addTarget(transmissionProcessList.back());
//        }
//    }
//
//    ModelTwo::Likelihood ModelTwo::value() {
//        return likelihood.value();
//    }
//
//    bool ModelTwo::isDirty() {
//        return likelihood.isDirty();
//    }
//
//
//    ModelTwo::State::State(std::map<std::string, LocusImpl *>& loci,
//                           std::vector<InfectionEvent *>& infections,
//                           std::map<InfectionEvent *, std::vector<InfectionEvent *>>& disallowedParents) : loci(loci), infections(infections), disallowedParents(disallowedParents) {
//        for (const auto &[locus_label, locus] : loci) {
//            alleleFrequencies.addLocus(locus);
//        }
//
//        infectionEventOrdering.addElements(infections);
//        observationFalsePositiveRate.initializeValue(.001);
//        observationFalseNegativeRate.initializeValue(.1);
//        geometricGenerationProb.initializeValue(.6);
//        lossProb.initializeValue(.9);
//        mutationProb.initializeValue(0.0);
//        mutationProb.initializeValue(1e-8);
//
//        meanCOI.initializeValue(3);
//    }
//}
//
