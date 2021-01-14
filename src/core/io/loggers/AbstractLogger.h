//
// Created by Maxwell Murphy on 1/12/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_ABSTRACTLOGGER_H
#define TRANSMISSION_NETWORKS_APP_ABSTRACTLOGGER_H

#include "AbstractOutput.h"

namespace transmission_nets::core::io {
    class AbstractLogger {
    public:

        template<typename Output>
        explicit AbstractLogger(std::unique_ptr<Output> output) : output_(std::move(output)) {}

        virtual std::string prepareValue() noexcept = 0;

        void log() {
            output_->write(prepareValue());
        }

    protected:
        std::unique_ptr<AbstractOutput> output_;
    };
}


#endif//TRANSMISSION_NETWORKS_APP_ABSTRACTLOGGER_H
