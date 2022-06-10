//
// Created by Maxwell Murphy on 8/18/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_CHAIN_H
#define TRANSMISSION_NETWORKS_APP_CHAIN_H

#include "core/io/loggers/FileOutput.h"
#include <vector>

namespace transmission_nets::core::mcmc {

    template<typename Model, typename SamplingScheduler, typename StateLogger>
    class Chain {
    public:
        void sample();
        void log();
        int totalSamples();

        [[nodiscard]] SamplingScheduler& getScheduler() const;
        [[nodiscard]] Model& getModel() const;
        [[nodiscard]] State& getState() const;
        [[nodiscard]] StateLogger& getStateLogger() const;

    private:
        Model model_;
        typename Model::State state_;
        SamplingScheduler scheduler_;
        StateLogger stateLogger_;
    };


    template<typename Model, typename SamplingScheduler, typename StateLogger>
    void Chain<Model, SamplingScheduler, StateLogger>::sample() {
        scheduler_.step();
    }


    template<typename Model, typename SamplingScheduler, typename StateLogger>
    void Chain<Model, SamplingScheduler, StateLogger>::log() {
        stateLogger_.log();
    }

    template<typename Model, typename SamplingScheduler, typename StateLogger>
    int Chain<Model, SamplingScheduler, StateLogger>::totalSamples() {
        return totalSamples_;
    }

    template<typename Model, typename SamplingScheduler, typename StateLogger>
    Model& Chain<Model, SamplingScheduler, StateLogger>::getModel() const {
        return model_;
    }

    template<typename Model, typename SamplingScheduler, typename StateLogger>
    State& Chain<Model, SamplingScheduler, StateLogger>::getState() const {
        return state_;
    }

    template<typename Model, typename SamplingScheduler, typename StateLogger>
    StateLogger& Chain<Model, SamplingScheduler, StateLogger>::getStateLogger() const {
        return stateLogger_;
    }

    template<typename Model, typename SamplingScheduler, typename StateLogger>
    SamplingScheduler& Chain<Model, SamplingScheduler, StateLogger>::getScheduler() const {
        return scheduler_;
    }

}// namespace transmission_nets::core::mcmc


#endif//TRANSMISSION_NETWORKS_APP_CHAIN_H