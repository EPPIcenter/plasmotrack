//
// Created by mmurphy on 10/29/21.
//

#include "StateLogger.h"

#include "core/io/loggers/CompressedFileOutput.h"
#include "core/io/loggers/ValueLogger.h"


namespace transmission_nets::impl::ModelNine {
    StateLogger::StateLogger(std::shared_ptr<State> state, fs::path rootPath, bool resetOutput) : state_(std::move(state)), rootPath_(std::move(rootPath)) {
        paramOutputFolder_ = rootPath_ / "parameters";
        auto epsPosFolder  = paramOutputFolder_ / "eps_pos";
        auto epsNegFolder  = paramOutputFolder_ / "eps_neg";
        auto infDurFolder  = paramOutputFolder_ / "infection_duration";
        auto freqDir       = paramOutputFolder_ / "allele_frequencies";


        if (!exists(paramOutputFolder_)) {
            create_directories(paramOutputFolder_);
        }

        if (!exists(epsPosFolder)) {
            create_directories(epsPosFolder);
        }

        if (!exists(epsNegFolder)) {
            create_directories(epsNegFolder);
        }

        if (!exists(infDurFolder)) {
            create_directories(infDurFolder);
        }

        if (!exists(freqDir)) {
            create_directories(freqDir);
        }

//        loggers_.push_back(new core::io::ValueLogger(state_->geometricGenerationProb, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "geo_gen_prob.csv", "geo_gen_prob", resetOutput)));
//        loggers_.push_back(new core::io::ValueLogger(state_->lossProb, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "loss_prob.csv", "loss_prob", resetOutput)));
//        loggers_.push_back(new core::io::ValueLogger(state_->mutationProb, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "mutation_prob.csv", "mutation_prob", resetOutput)));
        loggers_.push_back(new core::io::ValueLogger(state_->meanCOI, std::make_unique<core::io::CompressedFileOutput>(paramOutputFolder_ / "mean_coi.csv.gz", "mean_coi", resetOutput)));
        loggers_.push_back(new core::io::ValueLogger(state_->meanStrainsTransmitted, std::make_unique<core::io::CompressedFileOutput>(paramOutputFolder_ / "mean_strains_tx.csv.gz", "mean_strains_tx", resetOutput)));
        loggers_.push_back(new core::io::ValueLogger(state_->parentSetSizeProb, std::make_unique<core::io::CompressedFileOutput>(paramOutputFolder_ / "parent_set_size_prob.csv.gz", "parent_set_size_prob", resetOutput)));
//        loggers_.push_back(new core::io::ValueLogger(state_->symptomaticInfectionDurationShape, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "symptomatic_infection_duration_shape.csv", "infection_duration_shape", resetOutput)));
//        loggers_.push_back(new core::io::ValueLogger(state_->symptomaticInfectionDurationScale, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "symptomatic_infection_duration_scale.csv", "infection_duration_scale", resetOutput)));
//        loggers_.push_back(new core::io::ValueLogger(state_->asymptomaticInfectionDurationShape, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "asymptomatic_infection_duration_shape.csv", "infection_duration_shape", resetOutput)));
//        loggers_.push_back(new core::io::ValueLogger(state_->asymptomaticInfectionDurationScale, std::make_unique<core::io::FileOutput>(paramOutputFolder_ / "asymptomatic_infection_duration_scale.csv", "infection_duration_scale", resetOutput)));
        loggers_.push_back(new core::io::ValueLogger(state_->infectionEventOrdering, std::make_unique<core::io::CompressedFileOutput>(paramOutputFolder_ / "infection_order.csv.gz", "infection_order", resetOutput)));

        int i = 0;
        for (auto& infection : state_->infections) {
            auto inf_file_name = core::io::makePathValid(infection->id()) + ".csv.gz";
            loggers_.push_back(new core::io::ValueLogger(state_->expectedFalsePositives[i], std::make_unique<core::io::CompressedFileOutput>(epsPosFolder / inf_file_name, "eps_pos", resetOutput)));
            loggers_.push_back(new core::io::ValueLogger(state_->expectedFalseNegatives[i], std::make_unique<core::io::CompressedFileOutput>(epsNegFolder / inf_file_name, "eps_neg", resetOutput)));
            loggers_.push_back(new core::io::ValueLogger(infection->infectionDuration(), std::make_unique<core::io::CompressedFileOutput>(infDurFolder / inf_file_name, "duration", resetOutput)));
            i++;
        }

        for (const auto& [locus_label, locus] : state_->loci) {
            loggers_.push_back(new core::io::ValueLogger(state_->alleleFrequencies->alleleFrequencies(locus), std::make_unique<core::io::CompressedFileOutput>(freqDir / (locus_label + ".csv.gz"), "", resetOutput)));
        }

        for (const auto& infection : state_->infections) {
            auto inf_dir      = core::io::makePathValid(infection->id());
            auto full_inf_dir = paramOutputFolder_ / "genotypes" / inf_dir;
            if (!exists(full_inf_dir)) {
                create_directories(full_inf_dir);
            }
            for (const auto& [locus_label, locus] : state_->loci) {
                if (std::ranges::find(infection->loci(), locus) != infection->loci().end()) {
                    loggers_.push_back(new core::io::ValueLogger(infection->latentGenotype(locus), std::make_unique<core::io::CompressedFileOutput>(full_inf_dir / (locus_label + ".csv.gz"), "", resetOutput)));
                }
            }
        }

        for (const auto& infection : state_->latentParents) {
            auto inf_dir      = core::io::makePathValid(infection->id());
            auto full_inf_dir = paramOutputFolder_ / "latent_parents" / inf_dir;
            if (!exists(full_inf_dir)) {
                create_directories(full_inf_dir);
            }
            for (const auto& [locus_label, locus] : state_->loci) {
                if (std::ranges::find(infection->loci(), locus) != infection->loci().end()) {
                    loggers_.push_back(new core::io::ValueLogger(infection->latentGenotype(locus), std::make_unique<core::io::CompressedFileOutput>(full_inf_dir / (locus_label + ".csv.gz"), "", resetOutput)));
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
}// namespace transmission_nets::impl::ModelNine
