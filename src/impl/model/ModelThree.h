//
// Created by Maxwell Murphy on 6/18/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MODELTHREE_H
#define TRANSMISSION_NETWORKS_APP_MODELTHREE_H

#include "core/computation/Accumulator.h"

#include "core/distributions/ZTGeometric.h"
#include "core/distributions/ZTPoisson.h"
#include "core/distributions/ZTMultiplicativeBinomial.h"

#include "impl/state/ModelThreeState.h"

#include "model/observation_process/AlleleCounter.h"
#include "model/observation_process/ObservationProcessLikelihood.h"

#include "model/transmission_process/NetworkBasedTransmissionProcess.h"
#include "model/transmission_process/node_transmission_process/NoSuperInfectionMutation.h"

#include "model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h"


class ModelThree {
    static constexpr int MAX_COI = 10;
    static constexpr int MAX_PARENTS = 1;
    static constexpr int MAX_TRANSMISSIONS = 5;

    using GeneticsImpl = ModelThreeState::GeneticsImpl;
    using InfectionEvent = ModelThreeState::InfectionEvent;
    using AlleleFrequencyContainerImpl = ModelThreeState::AlleleFrequencyContainerImpl;

    using AlleleCounterImpl = AlleleCounter<GeneticsImpl>;
    using AlleleCounterAccumulator = Accumulator<AlleleCounterImpl, AlleleCounts>;

    using InterTransmissionProbImpl = ZTGeometric<MAX_TRANSMISSIONS>;
    using NodeTransmissionImpl = NoSuperInfectionMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>;

    using COIProbabilityImpl = ZTPoisson<MAX_COI>;
    using SourceTransmissionImpl = MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainerImpl, InfectionEvent, MAX_COI>;

    using ParentSetImpl = ParentSet<InfectionEvent>;
    using TransmissionProcess = NetworkBasedTransmissionProcess<MAX_PARENTS, NodeTransmissionImpl, SourceTransmissionImpl, InfectionEvent>;

//    using BetaPrior = Prior<boost::math::beta_distribution<>, Parameter<double>, int, int>;
//    using GammaPrior = Prior<boost::math::gamma_distribution<>, Parameter<double>, double, double>;

public:
    explicit ModelThree(ModelThreeState& state);

    double value();

    bool isDirty();

    ModelThreeState& state;
    Accumulator<PartialLikelihood, double> likelihood;

    // Observation Process
    std::vector<AlleleCounterImpl *> alleleCounters{};
    AlleleCounterAccumulator alleleCountAccumulator;
    ObservationProcessLikelihood<AlleleCounterAccumulator>* observationProcessLikelihood{};

    // Node Transmission Process
    InterTransmissionProbImpl* intp{};
    NodeTransmissionImpl* nodeTransmissionProcess{};

    // Source Transmission Process
    COIProbabilityImpl* coiProb{};
    std::vector<SourceTransmissionImpl *> sourceTransmissionProcessList;

    // Transmission Process
    std::vector<TransmissionProcess *> transmissionProcessList{};

};

#endif //TRANSMISSION_NETWORKS_APP_MODELTHREE_H
