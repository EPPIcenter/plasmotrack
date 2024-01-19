////
//// Created by Maxwell Murphy on 6/18/20.
////
//
//#include "ModelThree.h"
//
//#include "core/distributions/pdfs/BetaLogPDF.h"
//#include "core/distributions/pdfs/GammaLogPDF.h"
//
//namespace transmission_nets::impl {
//
//    ModelThree::ModelThree(ModelThreeState& state) : state(state) {
//        intp = new InterTransmissionProbImpl(state.geometricGenerationProb);
//        nodeTransmissionProcess = new NodeTransmissionImpl(state.mutationProb, state.lossProb, *intp);
//        coiProb = new COIProbabilityImpl(state.meanCOI);
//
//
//        // Register Priors
//        likelihood.addTarget(new core::distributions::BetaLogPDF(state.observationFalsePositiveRate, 1, 10));
//        likelihood.addTarget(new core::distributions::BetaLogPDF(state.observationFalseNegativeRate, 1, 10));
//        likelihood.addTarget(new core::distributions::BetaLogPDF(state.mutationProb, 1, 50));
//        likelihood.addTarget(new core::distributions::BetaLogPDF(state.lossProb, 1, 1));
//        likelihood.addTarget(new core::distributions::GammaLogPDF(state.meanCOI, 1, 1));
//        likelihood.addTarget(new core::distributions::BetaLogPDF(state.geometricGenerationProb, 1, 1));
//
//        for (auto &infection : state.infections) {
//            for (auto &[locus, obsGenotype] : infection->observedGenotype()) {
//                alleleCounters.push_back(new AlleleCounterImpl(infection->latentGenotype(locus), obsGenotype));
//                alleleCountAccumulator.addTarget(alleleCounters.back());
//            }
//        }
//
//        observationProcessLikelihood = new model::observation_process::ObservationProcessLikelihoodv1(
//                alleleCountAccumulator,
//                state.observationFalsePositiveRate,
//                state.observationFalseNegativeRate
//        );
//        likelihood.addTarget(observationProcessLikelihood);
//
//        for (unsigned int j = 0; j < state.infections.size(); ++j) {
//            auto infection = state.infections[j];
//            auto parentSet = state.transmissionNetwork.parentSet(infection);
//
//            sourceTransmissionProcessList.push_back(new SourceTransmissionImpl(
//                    *coiProb,
//                    state.alleleFrequencies,
//                    *infection
//                                                    )
//            );
//
//            transmissionProcessList.push_back(new TransmissionProcess(
//                    *nodeTransmissionProcess,
//                    *sourceTransmissionProcessList.back(),
//                    *infection,
//                    *parentSet
//                                              )
//            );
//
//            likelihood.addTarget(transmissionProcessList.back());
//
//        }
//    }
//
//    double ModelThree::value() {
//        return likelihood.value();
//    }
//
//    bool ModelThree::isDirty() {
//        return likelihood.isDirty();
//    }
//
//}
