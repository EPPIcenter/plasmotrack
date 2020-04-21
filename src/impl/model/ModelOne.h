//
// Created by Maxwell Murphy on 4/19/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MODELONE_H
#define TRANSMISSION_NETWORKS_APP_MODELONE_H

#include <utility>

#include "core/datatypes/Alleles.h"

#include "core/computation/PartialLikelihood.h"
#include "core/computation/Accumulator.h"

#include "impl/state/ModelOneState.h"

#include "model/observation_process/AlleleCounter.h"
#include "model/observation_process/AlleleCounts.h"
#include "model/observation_process/ObservationProcessLikelihood.h"

#include "model/transmission_process/OrderBasedTransmissionProcess.h"
#include "model/transmission_process/node_transmission_process/ZTMultiplicativeBinomial.h"
#include "model/transmission_process/node_transmission_process/NoSuperInfection.h"
#include "model/transmission_process/node_transmission_process/GeometricGenerationProbability.h"
#include "model/transmission_process/source_transmission_process/GeometricCOIProbability.h"
#include "model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h"

class ModelOne {
    static constexpr int MAX_COI = 10;
    static constexpr int MAX_PARENTS = 1;
    static constexpr int MAX_TRANSMISSIONS = 5;

    using GeneticsImpl = ModelOneState::GeneticsImpl;
    using InfectionEvent = ModelOneState::InfectionEvent;
    using AlleleFrequencyContainer = ModelOneState::AlleleFrequencyContainer;

    using AlleleCounter = AlleleCounter<GeneticsImpl>;
    using AlleleCounterAccumulator = Accumulator<AlleleCounter, AlleleCounts>;

    using Ordering = Ordering<InfectionEvent>;

    using COITransitionProbImpl = ZTMultiplicativeBinomial<MAX_COI>;
    using InterTransmissionProbImpl = GeometricGenerationProbability<MAX_TRANSMISSIONS>;
    using NodeTransmissionImpl = NoSuperInfection<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>;

    using COIProbabilityImpl = GeometricCOIProbability<MAX_COI>;
    using SourceTransmissionImpl = MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEvent>;

    using ParentSetImpl = OrderDerivedParentSet<InfectionEvent>;
    using TransmissionProcess = OrderBasedTransmissionProcess<MAX_PARENTS, NodeTransmissionImpl, SourceTransmissionImpl, InfectionEvent>;

public:
    explicit ModelOne(ModelOneState state);

    double value();

    void init();

    ModelOneState state;
    Accumulator<PartialLikelihood, double> likelihood;

    // Observation Process
    std::vector<std::unique_ptr<AlleleCounter>> alleleCounters{};
    AlleleCounterAccumulator alleleCountAccumulator;
    ObservationProcessLikelihood<AlleleCounterAccumulator>* observationProcessLikelihood{};

    // Node Transmission Process
    COITransitionProbImpl* coitp;
    InterTransmissionProbImpl* intp;
    NodeTransmissionImpl* nodeTransmissionProcess;

    // Source Transmission Process
    COIProbabilityImpl* coiProb;
    std::vector<std::unique_ptr<SourceTransmissionImpl>> sourceTransmissionProcessList;

    // Transmission Process
    std::vector<ParentSetImpl *> parentSetList{};
    std::vector<TransmissionProcess *> transmissionProcessList{};

public:


};


#endif //TRANSMISSION_NETWORKS_APP_MODELONE_H
