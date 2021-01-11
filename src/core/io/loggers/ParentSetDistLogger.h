//
// Created by Maxwell Murphy on 10/18/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_PARENTSETDISTLOGGER_H
#define TRANSMISSION_NETWORKS_APP_PARENTSETDISTLOGGER_H

#include <fstream>

#include "core/io/serialize.h"
#include "core/io/loggers/AbstractLogger.h"

namespace transmission_nets::core::io {

    template<typename TransmissionProcessImpl>
    class ParentSetDistLogger : public AbstractLogger {

    public:
        ParentSetDistLogger(fs::path outputPath, TransmissionProcessImpl* target);

        void logValue() noexcept override;

    private:
        TransmissionProcessImpl* target_;
        bool initialized{false};
        int iter = 1;

    };

    template<typename TransmissionProcessImpl>
    ParentSetDistLogger<TransmissionProcessImpl>::ParentSetDistLogger(fs::path outputPath, TransmissionProcessImpl* target) : AbstractLogger(outputPath), target_(target) {
        if(!fs::exists(outputPath_.parent_path())) {
            fs::create_directory(outputPath_.parent_path());
        }

        if(fs::exists(outputPath_)) {
            initialized = true;
        }

        outputFile_.open(outputPath, std::ofstream::app);
        outputFile_ << "parent_set,Llik,iter\n";
    }

    template<typename TransmissionProcessImpl>
    void ParentSetDistLogger<TransmissionProcessImpl>::logValue() noexcept {
        std::string out;
        const auto dist = target_->calcParentSetDist();
        for (const auto& [llik, ps] : dist.parentSetLliks) {
            out += serialize(ps) + "," + std::to_string(std::exp(llik - dist.totalLlik)) + "," + std::to_string(iter) + "\n";
        }
        out += "{S}," + std::to_string(std::exp(dist.sourceLlik - dist.totalLlik)) + "," + std::to_string(iter) + "\n";
        outputFile_ << out;
        iter++;
    }

}

#endif//TRANSMISSION_NETWORKS_APP_PARENTSETDISTLOGGER_H
