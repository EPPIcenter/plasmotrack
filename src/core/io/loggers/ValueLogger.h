//
// Created by Maxwell Murphy on 5/26/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_VALUELOGGER_H
#define TRANSMISSION_NETWORKS_APP_VALUELOGGER_H

#include <fstream>
#include <utility>

#include "AbstractLogger.h"
#include "core/io/serialize.h"

namespace transmission_nets::core::io {

    template<typename T>
    class ValueLogger : public AbstractLogger {

    public:
        template<typename Output>
        ValueLogger(std::shared_ptr<T> target, std::unique_ptr<Output> output);
        std::string prepareValue() noexcept override;

    private:
        std::shared_ptr<T> target_;
    };

    template<typename T>
    template<typename Output>
    ValueLogger<T>::ValueLogger(std::shared_ptr<T> target, std::unique_ptr<Output> output) : AbstractLogger(std::move(output)), target_(std::move(target)) {}

    template<typename T>
    std::string ValueLogger<T>::prepareValue() noexcept {
        return serialize(target_->value());
    }

}// namespace transmission_nets::core::io


#endif//TRANSMISSION_NETWORKS_APP_VALUELOGGER_H
