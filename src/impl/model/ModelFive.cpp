//
// Created by Maxwell Murphy on 4/13/21.
//

#include "ModelFive.h"

#include "core/distributions/pdfs/BetaLogPDF.h"
#include "core/distributions/pdfs/GammaLogPDF.h"
#include "core/io/loggers/FileOutput.h"
#include "core/io/loggers/MultiValueLogger.h"
#include "core/io/loggers/ParentSetDistLogger.h"
#include "core/io/loggers/ValueLogger.h"


namespace transmission_nets::impl::ModelFive {

    Model::Model(std::map<std::string, LocusImpl *>& loci,
                 std::vector<InfectionEvent *>& infections,
                 std::map<InfectionEvent *, std::vector<InfectionEvent *>>& disallowedParents) : state(loci, infections, disallowedParents) {
        intp = new InterTransmissionProbImpl(state.geometricGenerationProb);
        nodeTransmissionProcess = new NodeTransmissionImpl(state.mutationProb, state.lossProb, *intp);
        coiProb = new COIProbabilityImpl(state.meanCOI);

        // Register Priors
        likelihood.addTarget(new core::distributions::BetaLogPDF(state.mutationProb, state.mutationProbPriorAlpha, state.mutationProbPriorBeta));
        likelihood.addTarget(new core::distributions::BetaLogPDF(state.lossProb, state.lossProbPriorAlpha, state.lossProbPriorBeta));
        likelihood.addTarget(new core::distributions::BetaLogPDF(state.mutationProb, state.mutationProbPriorAlpha, state.mutationProbPriorBeta));
        likelihood.addTarget(new core::distributions::GammaLogPDF(state.meanCOI, state.meanCOIPriorShape, state.meanCOIPriorScale));
        likelihood.addTarget(new core::distributions::BetaLogPDF(state.geometricGenerationProb, state.geometricGenerationProbPriorAlpha, state.geometricGenerationProbPriorBeta));
//        likelihood.addTarget(new core::distributions::GammaLogPDF(state.infectionDurationShape, state.infectionDurationShapePriorShape, state.infectionDurationShapePriorScale));
//        likelihood.addTarget(new core::distributions::GammaLogPDF(state.infectionDurationScale, state.infectionDurationScalePriorShape, state.infectionDurationScalePriorScale));
        for (auto& obs : state.observationFalsePositiveRates) {
            likelihood.addTarget(new core::distributions::BetaLogPDF(obs, state.obsFPRPriorAlpha, state.obsFPRPriorBeta));
        }
        for (auto& obs : state.observationFalseNegativeRates) {
            likelihood.addTarget(new core::distributions::BetaLogPDF(obs, state.obsFNRPriorAlpha, state.obsFNRPriorBeta));
        }

        int i = 0;
        for (auto &infection : state.infections) {
            likelihood.addTarget(new core::distributions::GammaLogPDF(infection->infectionDuration(), state.infectionDurationShape, state.infectionDurationScale));
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
            parentSetList.push_back(new ParentSetImpl(&(state.infectionEventOrdering), infection));
            if (state.disallowedParents.contains(infection)) {
                parentSetList.back()->addDisallowedParents(state.disallowedParents.at(infection));
            }
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

    Likelihood Model::value() {
        return likelihood.value();
    }

    bool Model::isDirty() {
        return likelihood.isDirty();
    }


    State::State(std::map<std::string, LocusImpl *>& loci,
                 std::vector<InfectionEvent *>& infections,
                 std::map<InfectionEvent *, std::vector<InfectionEvent *>>& disallowedParents) : loci(loci), infections(infections), disallowedParents(disallowedParents) {
        obsFPRPriorAlpha.initializeValue(10);
        obsFPRPriorBeta.initializeValue(990);

        obsFNRPriorAlpha.initializeValue(10);
        obsFNRPriorBeta.initializeValue(900);

        geometricGenerationProbPriorAlpha.initializeValue(1);
        geometricGenerationProbPriorBeta.initializeValue(1);

        lossProbPriorAlpha.initializeValue(10);
        lossProbPriorBeta.initializeValue(90);

        mutationProbPriorAlpha.initializeValue(1);
        mutationProbPriorBeta.initializeValue(99);

        meanCOIPriorShape.initializeValue(1);
        meanCOIPriorScale.initializeValue(5);

//        infectionDurationShapePriorShape.initializeValue(1);
//        infectionDurationShapePriorScale.initializeValue(1000);
//        infectionDurationScalePriorShape.initializeValue(1);
//        infectionDurationScalePriorScale.initializeValue(1000);

        for (const auto &[locus_label, locus] : loci) {
            alleleFrequencies.addLocus(locus);
        }

        infectionEventOrdering.addElements(infections);

        infectionDurationShape.initializeValue(100);
        infectionDurationScale.initializeValue(1);

        for (size_t _ = 0; _ < infections.size(); ++_) {
            observationFalsePositiveRates.emplace_back();
            observationFalsePositiveRates.back().initializeValue(.01);
            observationFalseNegativeRates.emplace_back();
            observationFalseNegativeRates.back().initializeValue(.2);
        }


        geometricGenerationProb.initializeValue(.9);
        lossProb.initializeValue(.1);
        mutationProb.initializeValue(.01);
        meanCOI.initializeValue(5);
    }

    ParameterLogger::ParameterLogger(Model &model, fs::path rootPath) : model(model), rootPath(rootPath) {
        paramOutputFolder = rootPath / "parameters";
        statOutputFolder = rootPath / "stats";

        if (!fs::exists(paramOutputFolder)) {
            fs::create_directories(paramOutputFolder);
        }

        if (!fs::exists(statOutputFolder)) {
            fs::create_directories(statOutputFolder);
        }

        if (!fs::exists(paramOutputFolder / "eps_pos")) {
            fs::create_directories(paramOutputFolder / "eps_pos");
        }

        if (!fs::exists(paramOutputFolder / "eps_neg")) {
            fs::create_directories(paramOutputFolder / "eps_neg");
        }

        loggers.push_back(new core::io::ValueLogger(model.state.geometricGenerationProb, std::make_unique<core::io::FileOutput>(paramOutputFolder / "geo_gen_prob.csv", "geo_gen_prob")));
        loggers.push_back(new core::io::ValueLogger(model.state.lossProb, std::make_unique<core::io::FileOutput>(paramOutputFolder / "loss_prob.csv", "loss_prob")));
        loggers.push_back(new core::io::ValueLogger(model.state.mutationProb, std::make_unique<core::io::FileOutput>(paramOutputFolder / "mutation_prob.csv", "mutation_prob")));
        loggers.push_back(new core::io::ValueLogger(model.state.meanCOI, std::make_unique<core::io::FileOutput>(paramOutputFolder / "mean_coi.csv", "mean_coi")));
        loggers.push_back(new core::io::ValueLogger(model.state.infectionDurationShape, std::make_unique<core::io::FileOutput>(paramOutputFolder / "infection_duration_shape.csv", "infection_duration_shape")));
        loggers.push_back(new core::io::ValueLogger(model.state.infectionDurationScale, std::make_unique<core::io::FileOutput>(paramOutputFolder / "infection_duration_scale.csv", "infection_duration_scale")));
        loggers.push_back(new core::io::ValueLogger(model.state.infectionEventOrdering, std::make_unique<core::io::FileOutput>(paramOutputFolder / "infection_order.csv", "infection_order")));
        loggers.push_back(new core::io::ValueLogger(model, std::make_unique<core::io::FileOutput>(paramOutputFolder / "likelihood.csv", "llik")));


        auto eps_pos_logger = new core::io::MultiValueLogger<core::parameters::Parameter<double>>(std::make_unique<core::io::FileOutput>(paramOutputFolder / "eps_pos.csv", "id, eps_pos, iter"));
        auto eps_neg_logger = new core::io::MultiValueLogger<core::parameters::Parameter<double>>(std::make_unique<core::io::FileOutput>(paramOutputFolder / "eps_neg.csv", "id, eps_neg, iter"));
        auto infection_duration_logger = new core::io::MultiValueLogger<core::parameters::Parameter<double>>(std::make_unique<core::io::FileOutput>(paramOutputFolder / "infection_duration.csv", "id, infection_duration, iter"));
        int i = 0;
        for (auto &infection : model.state.infections) {
            eps_pos_logger->addTarget(infection->id(), model.state.observationFalsePositiveRates[i]);
            eps_neg_logger->addTarget(infection->id(), model.state.observationFalseNegativeRates[i]);
            infection_duration_logger->addTarget(infection->id(), infection->infectionDuration());
            i++;
        }
        loggers.push_back(eps_pos_logger);
        loggers.push_back(eps_neg_logger);
        loggers.push_back(infection_duration_logger);

        for (const auto &[locus_label, locus] : model.state.loci) {
            loggers.push_back(new core::io::ValueLogger(model.state.alleleFrequencies.alleleFrequencies(locus), std::make_unique<core::io::FileOutput>(paramOutputFolder / (locus_label + "_frequencies.csv"))));
        }

//        for (const auto &infection : model.state.infections) {
//            for (const auto &[locus_label, locus] : model.state.loci) {
//                if (std::find(infection->loci().begin(), infection->loci().end(), locus) != infection->loci().end()) {
//                    loggers.push_back(new core::io::ValueLogger(infection->latentGenotype(locus), std::make_unique<core::io::FileOutput>(paramOutputFolder / "nodes" / (infection->id() + "_" + locus_label + ".csv"))));
//                }
//            }
//        }

        for (const auto &tp : model.transmissionProcessList) {
            loggers.push_back(new core::io::ParentSetDistLogger(*tp, std::make_unique<core::io::FileOutput>(statOutputFolder / (tp->child_.id() + "_ps.csv"), "parent_set,Llik,iter")));
        }
    }

    void ParameterLogger::logParameters() const {
        for (const auto& logger : loggers) {
            logger->log();
        }
    }
}