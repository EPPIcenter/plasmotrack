//
// Created by Maxwell Murphy on 10/18/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_PARENTSETDISTLOGGER_H
#define TRANSMISSION_NETWORKS_APP_PARENTSETDISTLOGGER_H

#include <fstream>
#include <utility>

#include "AbstractLogger.h"
#include "core/io/serialize.h"

namespace transmission_nets::core::io {

    template<typename TransmissionProcessImpl>
    class ParentSetDistLogger : public AbstractLogger {

    public:
        template<typename Output>
        ParentSetDistLogger(TransmissionProcessImpl& target, std::unique_ptr<Output> output);
        std::string prepareValue() noexcept override;

    private:
        TransmissionProcessImpl& target_;
        int iter = 1;
    };

    template<typename TransmissionProcessImpl>
    template<typename Output>
    ParentSetDistLogger<TransmissionProcessImpl>::ParentSetDistLogger(TransmissionProcessImpl& target, std::unique_ptr<Output> output) : AbstractLogger(std::move(output)), target_(target) {}

    template<typename TransmissionProcessImpl>
    std::string ParentSetDistLogger<TransmissionProcessImpl>::prepareValue() noexcept {
        std::string out;
        const auto dist = target_.calcParentSetDist();
        for (const auto& [llik, ps] : dist.parentSetLliks) {
            out += serialize(ps) + "," + std::to_string(std::exp(llik - dist.totalLlik)) + "," + std::to_string(iter) + "\n";
        }
        out += "{S}," + std::to_string(std::exp(dist.sourceLlik - dist.totalLlik)) + "," + std::to_string(iter);
        iter++;
        return out;
    }

}

#endif//TRANSMISSION_NETWORKS_APP_PARENTSETDISTLOGGER_H
