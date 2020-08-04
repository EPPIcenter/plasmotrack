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

class ModelTwo {
    static constexpr int MAX_COI = 8;
    static constexpr int MAX_PARENTS = 1;
    static constexpr int MAX_TRANSMISSIONS = 5;

    using GeneticsImpl = ModelTwoState::GeneticsImpl;
    using InfectionEvent = ModelTwoState::InfectionEvent;
    using AlleleFrequencyContainer = ModelTwoState::AlleleFrequencyContainer;

    using AlleleCounter = AlleleCounter<GeneticsImpl>;
    using AlleleCounterAccumulator = Accumulator<AlleleCounter, AlleleCounts>;

    using Ordering = Ordering<InfectionEvent>;

    using InterTransmissionProbImpl = ZTGeometric<MAX_TRANSMISSIONS>;
    using NodeTransmissionImpl = NoSuperInfectionMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>;

    using COIProbabilityImpl = ZTPoisson<MAX_COI>;
    using SourceTransmissionImpl = MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEvent, MAX_COI>;

    using ParentSetImpl = OrderDerivedParentSet<InfectionEvent>;
    using TransmissionProcess = OrderBasedTransmissionProcess<MAX_PARENTS, NodeTransmissionImpl, SourceTransmissionImpl, InfectionEvent>;


public:
    explicit ModelTwo(ModelTwoState& state);

    double value();

    bool isDirty();

    ModelTwoState& state;
    Accumulator<PartialLikelihood, double> likelihood;

    // Observation Process
    std::vector<AlleleCounter *> alleleCounters{};
    AlleleCounterAccumulator alleleCountAccumulator;
    ObservationProcessLikelihood<AlleleCounterAccumulator>* observationProcessLikelihood{};

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


#endif //TRANSMISSION_NETWORKS_APP_MODELTWO_H
