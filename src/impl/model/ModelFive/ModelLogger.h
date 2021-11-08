//
// Created by mmurphy on 10/29/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_MODELLOGGER_H
#define TRANSMISSION_NETWORKS_APP_MODELLOGGER_H

#include "config.h"
#include "Model.h"

#include "core/io/loggers/ParentSetDistLogger.h"
#include "core/io/loggers/ValueLogger.h"

namespace transmission_nets::impl::ModelFive {

    class ModelLogger {
        std::shared_ptr<Model> model_;
        fs::path rootPath_;
        fs::path statOutputFolder_;
        std::vector<core::io::AbstractLogger*> loggers_{};

    public:
        ModelLogger(std::shared_ptr<Model> model, fs::path rootPath);
        ModelLogger(Model& model, fs::path rootPath);
        void log() const;
        void finalize() const;
    };
}



#endif//TRANSMISSION_NETWORKS_APP_MODELLOGGER_H
