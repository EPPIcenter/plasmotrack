//
// Created by mmurphy on 10/29/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_MODEL_H
#define TRANSMISSION_NETWORKS_APP_MODEL_H

#include "State.h"
#include "config.h"

#include <memory>

namespace transmission_nets::impl::ModelEight {
    struct Model : core::computation::PartialLikelihood {
        explicit Model(std::shared_ptr<State> state);
        explicit Model(State& state);

        Likelihood value() override;
        std::string identifier() override;

        std::shared_ptr<State> state_;

        core::computation::Accumulator<core::computation::PartialLikelihood, Likelihood> likelihood;

        // Observation Process
        std::vector<std::shared_ptr<model::observation_process::ObservationProcessLikelihoodv2<GeneticsImpl>>> observationProcessLikelihoodList{};

        // Node Transmission Process
        std::shared_ptr<InterTransmissionProbImpl> intp;
        std::shared_ptr<NodeTransmissionImpl> nodeTransmissionProcess;

        // Source Transmission Process
        std::shared_ptr<COIProbabilityImpl> coiProb;
        std::vector<std::shared_ptr<SourceTransmissionImpl>> sourceTransmissionProcessList{};

        // Transmission Process
        std::vector<std::shared_ptr<ParentSetImpl>> parentSetList{};
        std::vector<std::shared_ptr<TransmissionProcess>> transmissionProcessList{};
    };
}// namespace transmission_nets::impl::ModelEight


#endif//TRANSMISSION_NETWORKS_APP_MODEL_H
