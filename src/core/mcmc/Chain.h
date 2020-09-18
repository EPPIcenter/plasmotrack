//
// Created by Maxwell Murphy on 8/18/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_CHAIN_H
#define TRANSMISSION_NETWORKS_APP_CHAIN_H

#include <vector>

namespace transmission_nets::core::mcmc {

    template<typename Model, typename SamplingScheduler, typename Logger>
    class Chain {
    public:

        void sample() {
            scheduler_.step();
            totalSamples_++;
        };

        void log() {
            for (const auto& logger : loggers_) {
                logger->logValue();
            }
        };

        int totalSamples() {
            return totalSamples_;
        }


    //    void dumpState() {};
    //
    //    void restoreState() {};

    private:

        Model& model_;
        SamplingScheduler& scheduler_;
        std::vector<Logger*>&  loggers_;
        int totalSamples_ = 0;

    };

}



#endif//TRANSMISSION_NETWORKS_APP_CHAIN_H