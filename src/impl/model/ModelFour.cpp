////
//// Created by Maxwell Murphy on 12/10/20.
////
//
//#include "ModelFour.h"
//
//#include "core/distributions/pdfs/BetaLogPDF.h"
//#include "core/distributions/pdfs/GammaLogPDF.h"
//#include "core/io/loggers/ValueLogger.h"
//#include "core/io/loggers/ParentSetDistLogger.h"
//#include "core/io/loggers/FileOutput.h"
//
//
//namespace transmission_nets::impl::ModelFour {
//
//    Model::Model(std::map<std::string, LocusImpl *>& loci,
//                       std::vector<InfectionEvent *>& infections,
//                       std::map<InfectionEvent *, std::vector<InfectionEvent *>>& allowedParents) : state(loci, infections, allowedParents) {
//        intp = new InterTransmissionProbImpl(state.geometricGenerationProb);
//        nodeTransmissionProcess = new NodeTransmissionImpl(state.mutationProb, state.lossProb, *intp);
//        coiProb = new COIProbabilityImpl(state.meanCOI);
//
//        // Register Priors
//        likelihood.addTarget(new core::distributions::BetaLogPDF(state.mutationProb, state.mutationProbPriorAlpha, state.mutationProbPriorBeta));
//        likelihood.addTarget(new core::distributions::BetaLogPDF(state.lossProb, state.lossProbPriorAlpha, state.lossProbPriorBeta));
//        likelihood.addTarget(new core::distributions::BetaLogPDF(state.mutationProb, state.mutationProbPriorAlpha, state.mutationProbPriorBeta));
//        likelihood.addTarget(new core::distributions::GammaLogPDF(state.meanCOI, state.meanCOIPriorShape, state.meanCOIPriorScale));
//        likelihood.addTarget(new core::distributions::BetaLogPDF(state.geometricGenerationProb, state.geometricGenerationProbPriorAlpha, state.geometricGenerationProbPriorBeta));
//        for (auto& obs : state.observationFalsePositiveRates) {
//            likelihood.addTarget(new core::distributions::BetaLogPDF(obs, state.obsFPRPriorAlpha, state.obsFPRPriorBeta));
//        }
//        for (auto& obs : state.observationFalseNegativeRates) {
//            likelihood.addTarget(new core::distributions::BetaLogPDF(obs, state.obsFNRPriorAlpha, state.obsFNRPriorBeta));
//        }
//
//
//        int i = 0;
//        for (auto &infection : state.infections) {
//            alleleCountAccumulators.push_back(new AlleleCounterAccumulator());
//            for (auto &[locus, obsGenotype] : infection->observedGenotype()) {
//                alleleCounters.push_back(new AlleleCounterImpl(infection->latentGenotype(locus), obsGenotype));
//                alleleCountAccumulators.back()->addTarget(alleleCounters.back());
//            }
//            observationProcessLikelihood = new model::observation_process::ObservationProcessLikelihoodv1(
//                    *(alleleCountAccumulators.back()),
//                    state.observationFalseNegativeRates[i],
//                    state.observationFalsePositiveRates[i]);
//            likelihood.addTarget(observationProcessLikelihood);
//            parentSetList.push_back(new ParentSetImpl(&(state.infectionEventOrdering), infection));
//            if (state.allowedParents.contains(infection)) {
//                parentSetList.back()->addAllowedParents(state.allowedParents.at(infection));
//            }
//            i++;
//        }
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
//    Likelihood Model::value() {
//        return likelihood.value();
//    }
//
//    bool Model::isDirty() {
//        return likelihood.isDirty();
//    }
//
//
//    State::State(std::map<std::string, LocusImpl *>& loci,
//                           std::vector<InfectionEvent *>& infections,
//                           std::map<InfectionEvent *, std::vector<InfectionEvent *>>& allowedParents) : loci(loci), infections(infections), allowedParents(allowedParents) {
//        obsFPRPriorAlpha.initializeValue(10);
//        obsFPRPriorBeta.initializeValue(90);
//
//        obsFNRPriorAlpha.initializeValue(1);
//        obsFNRPriorBeta.initializeValue(99);
//
//        geometricGenerationProbPriorAlpha.initializeValue(1);
//        geometricGenerationProbPriorBeta.initializeValue(1);
//
//        lossProbPriorAlpha.initializeValue(10);
//        lossProbPriorBeta.initializeValue(90);
//
//        mutationProbPriorAlpha.initializeValue(1);
//        mutationProbPriorBeta.initializeValue(99);
//
//        meanCOIPriorShape.initializeValue(1);
//        meanCOIPriorShape.initializeValue(5);
//
//        for (const auto &[locus_label, locus] : loci) {
//            alleleFrequencies.addLocus(locus);
//        }
//
//        infectionEventOrdering.addElements(infections);
//
//        for (size_t _ = 0; _ < infections.size(); ++_) {
//            observationFalsePositiveRates.emplace_back();
//            observationFalsePositiveRates.back().initializeValue(.01);
//            observationFalseNegativeRates.emplace_back();
//            observationFalseNegativeRates.back().initializeValue(.2);
//        }
//
////        observationFalsePositiveRate.initializeValue(.001);
////        observationFalseNegativeRate.initializeValue(.1);
//
//        geometricGenerationProb.initializeValue(.9);
//        lossProb.initializeValue(.1);
//        mutationProb.initializeValue(.01);
//        meanCOI.initializeValue(5);
//    }
//
//    StateLogger::StateLogger(Model &model, fs::path rootPath) : model(model), rootPath(rootPath) {
//        paramOutputFolder = rootPath / "parameters";
//        statOutputFolder = rootPath / "stats";
//
//        if (!fs::exists(paramOutputFolder)) {
//            fs::create_directories(paramOutputFolder);
//        }
//
//        if (!fs::exists(statOutputFolder)) {
//            fs::create_directories(statOutputFolder);
//        }
//
//        if (!fs::exists(paramOutputFolder / "eps_pos")) {
//            fs::create_directories(paramOutputFolder / "eps_pos");
//        }
//
//        if (!fs::exists(paramOutputFolder / "eps_neg")) {
//            fs::create_directories(paramOutputFolder / "eps_neg");
//        }
//
//        loggers.push_back(new core::io::ValueLogger(model.state.geometricGenerationProb, std::make_unique<core::io::FileOutput>(paramOutputFolder / "geo_gen_prob.csv")));
//        loggers.push_back(new core::io::ValueLogger(model.state.lossProb, std::make_unique<core::io::FileOutput>(paramOutputFolder / "loss_prob.csv")));
//        loggers.push_back(new core::io::ValueLogger(model.state.mutationProb, std::make_unique<core::io::FileOutput>(paramOutputFolder / "mutation_prob.csv")));
//        loggers.push_back(new core::io::ValueLogger(model.state.meanCOI, std::make_unique<core::io::FileOutput>(paramOutputFolder / "mean_coi.csv")));
//        loggers.push_back(new core::io::ValueLogger(model.state.infectionEventOrdering, std::make_unique<core::io::FileOutput>(paramOutputFolder / "infection_order.csv")));
//        loggers.push_back(new core::io::ValueLogger(model, std::make_unique<core::io::FileOutput>(paramOutputFolder / "likelihood.csv")));
//
//        int i = 0;
//        for (auto &infection : model.state.infections) {
//            loggers.push_back(new core::io::ValueLogger(model.state.observationFalsePositiveRates[i], std::make_unique<core::io::FileOutput>(paramOutputFolder / "eps_pos" / (infection->id() + ".csv"))));
//            loggers.push_back(new core::io::ValueLogger(model.state.observationFalseNegativeRates[i], std::make_unique<core::io::FileOutput>(paramOutputFolder / "eps_neg" / (infection->id() + ".csv"))));
//            i++;
//        }
//
//        for (const auto &[locus_label, locus] : model.state.loci) {
//            loggers.push_back(new core::io::ValueLogger(model.state.alleleFrequencies.alleleFrequencies(locus), std::make_unique<core::io::FileOutput>(paramOutputFolder / (locus_label + "_frequencies.csv"))));
//        }
//
//        for (const auto &infection : model.state.infections) {
//            for (const auto &[locus_label, locus] : model.state.loci) {
//                if (std::find(infection->loci().begin(), infection->loci().end(), locus) != infection->loci().end()) {
//                    loggers.push_back(new core::io::ValueLogger(infection->latentGenotype(locus), std::make_unique<core::io::FileOutput>(paramOutputFolder / "nodes" / (infection->id() + "_" + locus_label + ".csv"))));
//                }
//            }
//        }
//
//        for (const auto &tp : model.transmissionProcessList) {
//            loggers.push_back(new core::io::ParentSetDistLogger(*tp, std::make_unique<core::io::FileOutput>(statOutputFolder / (tp->child_.id() + "_ps.csv"), "parent_set,Llik,iter")));
//        }
//    }
//
//    void StateLogger::logParameters() const {
//        for (const auto& logger : loggers) {
//            logger->log();
//        }
//    }
//}