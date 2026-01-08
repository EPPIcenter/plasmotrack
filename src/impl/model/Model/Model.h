//
// Created by mmurphy on 10/29/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_MODEL_H
#define TRANSMISSION_NETWORKS_APP_MODEL_H

#include "State.h"
#include "config.h"
#include "core/computation/Accumulator.h"

#include <memory>
#include <shared_mutex>

namespace transmission_nets::impl::Model {
    struct Model : core::computation::PartialLikelihood {
        using ModelLikelihood = core::computation::Accumulator<PartialLikelihood, Likelihood>;
        explicit Model(std::shared_ptr<State> state, double temperature = 1.0);
        explicit Model(State& state, double temperature = 1.0);

        Likelihood value() override;
        Likelihood valueThreadSafe();
        std::string identifier() override;
        [[nodiscard]] double getTemperature() const;
        void setTemperature(double t);

        [[nodiscard]] double getPrior();
        double getLikelihood();

        std::shared_ptr<State> state_;

        ModelLikelihood likelihood;
        ModelLikelihood prior;
        double temperature;

        // Observation Process
        std::vector<std::shared_ptr<model::observation_process::ObservationProcessLikelihoodv2<GeneticsImpl>>> observationProcessLikelihoodList{};

        // Parent Set Size Likelihood
        std::shared_ptr<ParentSetSizeLikelihoodImpl> parentSetSizeLikelihood;

        // Node Transmission Process
//        std::shared_ptr<InterTransmissionProbImpl> intp;
        std::shared_ptr<NodeTransmissionImpl> nodeTransmissionProcess;

        // Source Transmission Process
        std::shared_ptr<COIProbabilityImpl> coiProb;
        std::vector<std::shared_ptr<SourceTransmissionImpl>> sourceTransmissionProcessList{};

        // Transmission Process
//        std::map<std::shared_ptr<InfectionEvent>, std::shared_ptr<ParentSetImpl>> parentSetList{};
        std::vector<std::shared_ptr<TransmissionProcess>> transmissionProcessList{};
        std::shared_mutex valueMutex;

    };
}// namespace transmission_nets::impl::Model


#endif//TRANSMISSION_NETWORKS_APP_MODEL_H
