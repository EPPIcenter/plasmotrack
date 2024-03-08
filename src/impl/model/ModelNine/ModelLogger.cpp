//
// Created by mmurphy on 10/29/21.
//

#include "ModelLogger.h"


namespace transmission_nets::impl::ModelNine {
    ModelLogger::ModelLogger(std::shared_ptr<Model> model, fs::path rootPath) : model_(std::move(model)), rootPath_(std::move(rootPath)) {
        const auto statOutputFolder_ = rootPath_ / "stats";
        const auto parentSetFolder_  = rootPath_ / "parent_sets";
        if (!exists(statOutputFolder_)) {
            create_directories(statOutputFolder_);
        }

        if (!exists(parentSetFolder_)) {
            create_directories(parentSetFolder_);
        }

        loggers_.push_back(new core::io::ValueLogger(model_, std::make_unique<core::io::CompressedFileOutput>(statOutputFolder_ / "likelihood.csv.gz", "llik")));

        for (const auto& tp : model_->transmissionProcessList) {
            auto path = parentSetFolder_ / (core::io::makePathValid(tp->child_->id() + "_ps.csv.gz"));
            loggers_.push_back(new core::io::ParentSetDistLogger(tp, std::make_unique<core::io::CompressedFileOutput>(path, "parent_set,prob,iter"), path, false));
        }
    }

    void ModelLogger::log() const {
        for (const auto& logger : loggers_) {
            logger->log();
        }
    }

    void ModelLogger::finalize() const {
        for (const auto& logger : loggers_) {
            logger->finalize();
        }
    }
}// namespace transmission_nets::impl::ModelNine