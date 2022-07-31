//
// Created by mmurphy on 10/29/21.
//

#include "ModelLogger.h"


namespace transmission_nets::impl::ModelEight {
    ModelLogger::ModelLogger(std::shared_ptr<Model> model, fs::path rootPath) : model_(std::move(model)), rootPath_(std::move(rootPath)) {
        auto statOutputFolder_ = rootPath_ / "stats";
        auto parentSetFolder_  = rootPath_ / "parent_sets";
        if (!fs::exists(statOutputFolder_)) {
            fs::create_directories(statOutputFolder_);
        }

        if (!fs::exists(parentSetFolder_)) {
            fs::create_directories(parentSetFolder_);
        }

        loggers_.push_back(new core::io::ValueLogger(model_, std::make_unique<core::io::FileOutput>(statOutputFolder_ / "likelihood.csv", "llik")));

        for (const auto& tp : model_->transmissionProcessList) {
            auto path = parentSetFolder_ / (core::io::makePathValid(tp->child_->id() + "_ps.csv"));
            loggers_.push_back(new core::io::ParentSetDistLogger(tp, std::make_unique<core::io::FileOutput>(path, "parent_set,prob,iter"), path, false));
        }
    }

    ModelLogger::ModelLogger(Model& model, fs::path rootPath) : ModelLogger(std::make_shared<Model>(model), std::move(rootPath)) {}

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
}// namespace transmission_nets::impl::ModelEight