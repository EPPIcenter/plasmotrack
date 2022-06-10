//
// Created by mmurphy on 10/29/21.
//

#include "StateLogger.h"


namespace transmission_nets::impl::ModelSeven {
    StateLogger::StateLogger(std::shared_ptr<State> state, fs::path rootPath, bool resetOutput) : state_(std::move(state)), rootPath_(std::move(rootPath)) {
        paramOutputFolder_ = rootPath_ / "parameters";
        auto epsPosFolder  = paramOutputFolder_ / "eps_pos";
        auto epsNegFolder  = paramOutputFolder_ / "eps_neg";
        auto infDurFolder  = paramOutputFolder_ / "infection_duration";
        auto freqDir       = paramOutputFolder_ / "allele_frequencies";

        if (!fs::exists(paramOutputFolder_)) {
            fs::create_directories(paramOutputFolder_);
        }

        if (!fs::exists(epsPosFolder)) {
            fs::create_directories(epsPosFolder);
        }

        if (!fs::exists(epsNegFolder)) {
            fs::create_directories(epsNegFolder);
        }

        if (!fs::exists(infDurFolder)) {
            fs::create_directories(infDurFolder);
        }

        if (!fs::exists(freqDir)) {
            fs::create_directories(freqDir);
        }

        loggers_.push_back(new core::io::ValueLogger(state_->geometricGenerationProb, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "geo_gen_prob.csv", "geo_gen_prob", resetOutput)));
        loggers_.push_back(new core::io::ValueLogger(state_->lossProb, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "loss_prob.csv", "loss_prob", resetOutput)));
        loggers_.push_back(new core::io::ValueLogger(state_->mutationProb, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "mutation_prob.csv", "mutation_prob", resetOutput)));
        loggers_.push_back(new core::io::ValueLogger(state_->meanCOI, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "mean_coi.csv", "mean_coi", resetOutput)));
        loggers_.push_back(new core::io::ValueLogger(state_->infectionDurationShape, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "infection_duration_shape.csv", "infection_duration_shape", resetOutput)));
        loggers_.push_back(new core::io::ValueLogger(state_->infectionDurationScale, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "infection_duration_scale.csv", "infection_duration_scale", resetOutput)));
        loggers_.push_back(new core::io::ValueLogger(state_->infectionEventOrdering, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "infection_order.csv", "infection_order", resetOutput)));

        int i = 0;
        for (auto& infection : state_->infections) {
            auto inf_file_name = core::io::makePathValid(infection->id()) + ".csv";
            loggers_.push_back(new core::io::ValueLogger(state_->expectedFalsePositives[i], std::make_unique<core::io::FileOutput>(epsPosFolder / inf_file_name, "eps_pos", resetOutput)));
            loggers_.push_back(new core::io::ValueLogger(state_->expectedFalseNegatives[i], std::make_unique<core::io::FileOutput>(epsNegFolder / inf_file_name, "eps_neg", resetOutput)));
            loggers_.push_back(new core::io::ValueLogger(infection->infectionDuration(), std::make_unique<core::io::FileOutput>(infDurFolder / inf_file_name, "duration", resetOutput)));
            i++;
        }

        for (const auto& [locus_label, locus] : state_->loci) {
            loggers_.push_back(new core::io::ValueLogger(state_->alleleFrequencies->alleleFrequencies(locus), std::make_unique<core::io::FileOutput>(freqDir / (locus_label + ".csv"), "", resetOutput)));
        }

        for (const auto& infection : state_->infections) {
            auto inf_dir      = core::io::makePathValid(infection->id());
            auto full_inf_dir = paramOutputFolder_ / "genotypes" / inf_dir;
            if (!fs::exists(full_inf_dir)) {
                fs::create_directories(full_inf_dir);
            }
            for (const auto& [locus_label, locus] : state_->loci) {
                if (std::find(infection->loci().begin(), infection->loci().end(), locus) != infection->loci().end()) {
                    loggers_.push_back(new core::io::ValueLogger(infection->latentGenotype(locus), std::make_unique<core::io::FileOutput>(full_inf_dir / (locus_label + ".csv"), "", resetOutput)));
                }
            }
        }
    }


    void StateLogger::log() const {
        for (const auto& logger : loggers_) {
            logger->log();
        }
    }

    void StateLogger::finalize() const {
        for (const auto& logger : loggers_) {
            logger->finalize();
        }
    }
}// namespace transmission_nets::impl::ModelSeven
