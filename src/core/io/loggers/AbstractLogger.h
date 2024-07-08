//
// Created by Maxwell Murphy on 1/12/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_ABSTRACTLOGGER_H
#define TRANSMISSION_NETWORKS_APP_ABSTRACTLOGGER_H

#include "AbstractOutput.h"

#include <memory>

namespace transmission_nets::core::io {
    class AbstractLogger {
    public:
        virtual ~AbstractLogger() = default;
        template<typename Output>
        explicit AbstractLogger(std::unique_ptr<Output> output) : output_(std::move(output)) {}

        virtual std::string prepareValue() noexcept = 0;

        void log() {
            output_->write(prepareValue());
        }

        void finalize() {
            output_->finalize();
        }

    protected:
        std::unique_ptr<AbstractOutput> output_{};
    };
}// namespace transmission_nets::core::io


#endif//TRANSMISSION_NETWORKS_APP_ABSTRACTLOGGER_H
