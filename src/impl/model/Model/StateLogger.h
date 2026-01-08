//
// Created by mmurphy on 10/29/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_STATELOGGER_H
#define TRANSMISSION_NETWORKS_APP_STATELOGGER_H

#include "State.h"
#include "config.h"
#include "core/io/loggers/AbstractLogger.h"
#include "core/io/utils.h"

#include <memory>
#include <regex>
#include <utility>

namespace transmission_nets::impl::Model {
    class StateLogger {
        std::shared_ptr<State> state_;
        fs::path rootPath_;
        fs::path paramOutputFolder_;
        std::vector<core::io::AbstractLogger*> loggers_{};

    public:
        StateLogger(std::shared_ptr<State> state, fs::path rootPath, bool resetOutput = false);
        //        StateLogger(State& state, fs::path rootPath);
        void log() const;
        void finalize() const;
    };


}// namespace transmission_nets::impl::Model


#endif//TRANSMISSION_NETWORKS_APP_STATELOGGER_H
