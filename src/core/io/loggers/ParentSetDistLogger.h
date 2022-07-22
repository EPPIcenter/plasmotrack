//
// Created by Maxwell Murphy on 10/18/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_PARENTSETDISTLOGGER_H
#define TRANSMISSION_NETWORKS_APP_PARENTSETDISTLOGGER_H


#include "AbstractLogger.h"
#include "core/io/serialize.h"

#include <fmt/core.h>

#include <fstream>
#include <utility>


namespace transmission_nets::core::io {

    template<typename TransmissionProcessImpl>
    class ParentSetDistLogger : public AbstractLogger {

    public:
        template<typename Output>
        ParentSetDistLogger(std::shared_ptr<TransmissionProcessImpl> target, std::unique_ptr<Output> output, std::string  uid, bool include_source = true);
        std::string prepareValue() noexcept override;

    private:
        static std::map<std::string, int>& getIterCounters();
        void incrementIterCounter();
        int getIterCounter();

        std::shared_ptr<TransmissionProcessImpl> target_;
        std::string uid_;
        bool include_source_;
    };

    template<typename TransmissionProcessImpl>
    template<typename Output>
    ParentSetDistLogger<TransmissionProcessImpl>::ParentSetDistLogger(std::shared_ptr<TransmissionProcessImpl> target, std::unique_ptr<Output> output, std::string  uid, bool include_source) : AbstractLogger(std::move(output)), target_(target), uid_(std::move(uid)), include_source_(include_source){}

    template<typename TransmissionProcessImpl>
    std::string ParentSetDistLogger<TransmissionProcessImpl>::prepareValue() noexcept {
        std::string out;
        int iter_       = getIterCounter();
        const auto dist = target_->calcParentSetDist();
        for (const auto& [llik, ps] : dist.parentSetLliks) {
            out += fmt::format("{},{},{}\n", serialize(ps), std::exp(llik - dist.totalLlik), iter_);
        }
        if (include_source_) {
            out += fmt::format("{{S}},{},{}", std::exp(dist.sourceLlik - dist.totalLlik), iter_);
        }
        incrementIterCounter();
        return out;
    }

    template<typename TransmissionProcessImpl>
    std::map<std::string, int>& ParentSetDistLogger<TransmissionProcessImpl>::getIterCounters() {
        static std::map<std::string, int> iterCounters_{};
        return iterCounters_;
    }

    template<typename TransmissionProcessImpl>
    void ParentSetDistLogger<TransmissionProcessImpl>::incrementIterCounter() {
        auto& counters = getIterCounters();
        counters[uid_] = counters[uid_] + 1;
    }

    template<typename TransmissionProcessImpl>
    int ParentSetDistLogger<TransmissionProcessImpl>::getIterCounter() {
        return getIterCounters()[uid_];
    }


}// namespace transmission_nets::core::io

#endif//TRANSMISSION_NETWORKS_APP_PARENTSETDISTLOGGER_H
