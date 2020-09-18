//
// Created by Maxwell Murphy on 6/5/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MODELTWO_H
#define TRANSMISSION_NETWORKS_APP_MODELTWO_H

#include <utility>

#include "core/computation/PartialLikelihood.h"
#include "core/computation/Accumulator.h"
#include "core/computation/Transformers/LogTransformer.h"

#include "core/datatypes/Alleles.h"

#include "core/distributions/ZTGeometric.h"
#include "core/distributions/ZTPoisson.h"
#include "core/distributions/ZTMultiplicativeBinomial.h"

#include "impl/state/ModelTwoState.h"

#include "model/observation_process/AlleleCounter.h"
#include "model/observation_process/AlleleCounts.h"
#include "model/observation_process/ObservationProcessLikelihood.h"

#include "model/transmission_process/OrderBasedTransmissionProcess.h"
#include "model/transmission_process/node_transmission_process/NoSuperInfectionMutation.h"

#include "model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h"

namespace transmission_nets::impl {

    class ModelTwo {
        static constexpr int MAX_COI = 8;
        static constexpr int MAX_PARENTS = 1;
        static constexpr int MAX_TRANSMISSIONS = 5;

        using GeneticsImpl = ModelTwoState::GeneticsImpl;
        using InfectionEvent = ModelTwoState::InfectionEvent;
        using AlleleFrequencyContainerImpl = ModelTwoState::AlleleFrequencyContainerImpl;

        using AlleleCounterImpl = model::observation_process::AlleleCounter<GeneticsImpl>;
        using AlleleCounterAccumulator = core::computation::Accumulator<AlleleCounterImpl, model::observation_process::AlleleCounts>;

//    using OrderingImpl = Ordering<InfectionEvent>;

        using InterTransmissionProbImpl = core::distributions::ZTGeometric<MAX_TRANSMISSIONS>;
        using NodeTransmissionImpl = model::transmission_process::NoSuperInfectionMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>;

        using COIProbabilityImpl = core::distributions::ZTPoisson<MAX_COI>;
        using SourceTransmissionImpl = model::transmission_process::MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainerImpl, InfectionEvent, MAX_COI>;

        using ParentSetImpl = core::computation::OrderDerivedParentSet<InfectionEvent>;
        using TransmissionProcess = model::transmission_process::OrderBasedTransmissionProcess<MAX_PARENTS, NodeTransmissionImpl, SourceTransmissionImpl, InfectionEvent>;


    public:
        explicit ModelTwo(ModelTwoState& state);

        double value();

        bool isDirty();

        ModelTwoState& state;
        core::computation::Accumulator<core::computation::PartialLikelihood, double> likelihood;

        // Observation Process
        std::vector<AlleleCounterImpl *> alleleCounters{};
        AlleleCounterAccumulator alleleCountAccumulator;
        model::observation_process::ObservationProcessLikelihood<AlleleCounterAccumulator>* observationProcessLikelihood{};

        // Node Transmission Process
        InterTransmissionProbImpl* intp{};
        NodeTransmissionImpl* nodeTransmissionProcess{};

        // Source Transmission Process
        COIProbabilityImpl* coiProb{};
        std::vector<SourceTransmissionImpl *> sourceTransmissionProcessList;

        // Transmission Process
        std::vector<ParentSetImpl *> parentSetList{};
        std::vector<TransmissionProcess *> transmissionProcessList{};

    };

}


#endif //TRANSMISSION_NETWORKS_APP_MODELTWO_H
