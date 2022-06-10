//
// Created by Maxwell Murphy on 6/22/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_LAMBDALOGGER_H
#define TRANSMISSION_NETWORKS_APP_LAMBDALOGGER_H

#include <fstream>
#include <utility>

#include "AbstractLogger.h"

namespace transmission_nets::core::io {

    template<typename Loggable>
    class LambdaLogger : public AbstractLogger {

    public:
        template<typename Output>
        LambdaLogger(Loggable f, std::unique_ptr<Output> output);
        std::string prepareValue() noexcept override;

    private:
        Loggable f_;
    };

    template<typename Loggable>
    template<typename Output>
    LambdaLogger<Loggable>::LambdaLogger(Loggable f, std::unique_ptr<Output> output) : AbstractLogger(std::move(output)), f_(f) {}

    template<typename Loggable>
    std::string LambdaLogger<Loggable>::prepareValue() noexcept {
        return f_();
    }

}// namespace transmission_nets::core::io

#endif//TRANSMISSION_NETWORKS_APP_LAMBDALOGGER_H
