//
// Created by Maxwell Murphy on 4/13/21.
//

//#include "ModelFive.h"
//
//#include <utility>
//
//#include "core/distributions/pdfs/BetaLogPDF.h"
//#include "core/distributions/pdfs/GammaLogPDF.h"
//#include "core/io/loggers/FileOutput.h"
//#include "core/io/loggers/MultiValueLogger.h"
//#include "core/io/loggers/ParentSetDistLogger.h"
//#include "core/io/loggers/ValueLogger.h"
//
//
//namespace transmission_nets::impl::ModelFive {
//
//
//
//
//
//    StateLogger::StateLogger(State &state, fs::path rootPath) : state_(state), rootPath_(std::move(rootPath)) {
//        paramOutputFolder_ = rootPath_ / "parameters";
//
//        if (!fs::exists(paramOutputFolder_)) {
//            fs::create_directories(paramOutputFolder_);
//        }
//
//
//        if (!fs::exists(paramOutputFolder_ / "eps_pos")) {
//            fs::create_directories(paramOutputFolder_ / "eps_pos");
//        }
//
//        if (!fs::exists(paramOutputFolder_ / "eps_neg")) {
//            fs::create_directories(paramOutputFolder_ / "eps_neg");
//        }
//
//        loggers_.push_back(new core::io::ValueLogger(state_.geometricGenerationProb, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "geo_gen_prob.csv", "geo_gen_prob")));
//        loggers_.push_back(new core::io::ValueLogger(state_.lossProb, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "loss_prob.csv", "loss_prob")));
//        loggers_.push_back(new core::io::ValueLogger(state_.mutationProb, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "mutation_prob.csv", "mutation_prob")));
//        loggers_.push_back(new core::io::ValueLogger(state_.meanCOI, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "mean_coi.csv", "mean_coi")));
//        loggers_.push_back(new core::io::ValueLogger(state_.infectionDurationShape, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "infection_duration_shape.csv", "infection_duration_shape")));
//        loggers_.push_back(new core::io::ValueLogger(state_.infectionDurationScale, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "infection_duration_scale.csv", "infection_duration_scale")));
//        loggers_.push_back(new core::io::ValueLogger(state_.infectionEventOrdering, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "infection_order.csv", "infection_order")));
//
//        auto eps_pos_logger = new core::io::MultiValueLogger<core::parameters::Parameter<float>>(std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "eps_pos.csv", "id, eps_pos, iter"));
//        auto eps_neg_logger = new core::io::MultiValueLogger<core::parameters::Parameter<float>>(std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "eps_neg.csv", "id, eps_neg, iter"));
//        auto infection_duration_logger = new core::io::MultiValueLogger<core::parameters::Parameter<float>>(std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "infection_duration.csv", "id, infection_duration, iter"));
//        int i = 0;
//        for (auto &infection : state_.infections) {
//            eps_pos_logger->addTarget(infection->id(), state_.observationFalsePositiveRates[i]);
//            eps_neg_logger->addTarget(infection->id(), state_.observationFalseNegativeRates[i]);
//            infection_duration_logger->addTarget(infection->id(), infection->infectionDuration());
//            i++;
//        }
//        loggers_.push_back(eps_pos_logger);
//        loggers_.push_back(eps_neg_logger);
//        loggers_.push_back(infection_duration_logger);
//
//        for (const auto &[locus_label, locus] : state_.loci) {
//            loggers_.push_back(new core::io::ValueLogger(state_.alleleFrequencies.alleleFrequencies(locus), std::make_unique<core::io::FileOutput>(paramOutputFolder_ / (locus_label + "_frequencies.csv"))));
//        }
//
//        for (const auto &infection : state_.infections) {
//            for (const auto &[locus_label, locus] : state_.loci) {
//                if (std::find(infection->loci().begin(), infection->loci().end(), locus) != infection->loci().end()) {
//                    loggers_.push_back(new core::io::ValueLogger(infection->latentGenotype(locus), std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "nodes" / (infection->id() + "_" + locus_label + ".csv"))));
//                }
//            }
//        }
//
//
//    }
//
//    void StateLogger::log() const {
//        for (const auto& logger : loggers_) {
//            logger->log();
//        }
//    }
//
//    ModelLogger::ModelLogger(Model &model, fs::path rootPath) : model_(model), rootPath_(std::move(rootPath)) {
//        statOutputFolder_ = rootPath_ / "stats";
//        if (!fs::exists(statOutputFolder_)) {
//            fs::create_directories(statOutputFolder_);
//        }
//
//        loggers_.push_back(new core::io::ValueLogger(model_, std::make_unique<core::io::FileOutput>(statOutputFolder_ / "likelihood.csv", "llik")));
//
//        for (const auto &tp : model_.transmissionProcessList) {
//            loggers_.push_back(new core::io::ParentSetDistLogger(*tp, std::make_unique<core::io::FileOutput>(statOutputFolder_ / (tp->child_.id() + "_ps.csv"), "parent_set,Llik,iter")));
//        }
//    }
//
//    void ModelLogger::log() const {
//        for (const auto& logger : loggers_) {
//            logger->log();
//        }
//    }
//
//
//}