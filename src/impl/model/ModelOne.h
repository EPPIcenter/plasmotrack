//
// Created by Maxwell Murphy on 4/19/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MODELONE_H
#define TRANSMISSION_NETWORKS_APP_MODELONE_H

#include <utility>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/gamma.hpp>

#include "core/computation/PartialLikelihood.h"
#include "core/computation/Accumulator.h"
#include "core/computation/Transformers/LogTransformer.h"

#include "core/datatypes/Alleles.h"

#include "core/distributions/ZTGeometric.h"
#include "core/distributions/ZTMultiplicativeBinomial.h"

#include "core/priors/Prior.h"

#include "impl/state/ModelOneState.h"

#include "model/observation_process/AlleleCounter.h"
#include "model/observation_process/AlleleCounts.h"
#include "model/observation_process/ObservationProcessLikelihood.h"

#include "model/transmission_process/OrderBasedTransmissionProcess.h"
#include "model/transmission_process/node_transmission_process/NoSuperInfectionNoMutation.h"

#include "model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h"

class ModelOne {
    static constexpr int MAX_COI = 32;
    static constexpr int MAX_PARENTS = 1;
    static constexpr int MAX_TRANSMISSIONS = 5;

    using GeneticsImpl = ModelOneState::GeneticsImpl;
    using InfectionEvent = ModelOneState::InfectionEvent;
    using AlleleFrequencyContainer = ModelOneState::AlleleFrequencyContainer;

    using AlleleCounter = AlleleCounter<GeneticsImpl>;
    using AlleleCounterAccumulator = Accumulator<AlleleCounter, AlleleCounts>;

    using Ordering = Ordering<InfectionEvent>;

    using COITransitionProbImpl = ZTMultiplicativeBinomial<MAX_COI>;
    using InterTransmissionProbImpl = ZTGeometric<MAX_TRANSMISSIONS>;
    using NodeTransmissionImpl = NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>;

    using COIProbabilityImpl = ZTGeometric<MAX_COI>;
    using SourceTransmissionImpl = MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainer, InfectionEvent, MAX_COI>;

    using ParentSetImpl = OrderDerivedParentSet<InfectionEvent>;
    using TransmissionProcess = OrderBasedTransmissionProcess<MAX_PARENTS, NodeTransmissionImpl, SourceTransmissionImpl, InfectionEvent>;

    using BetaPrior = Prior<boost::math::beta_distribution<>, Parameter<double>, int, int>;
    using GammaPrior = Prior<boost::math::gamma_distribution<>, Parameter<double>, double, double>;

public:
    explicit ModelOne(ModelOneState& state);

    double value();

    bool isDirty();

    ModelOneState& state;
    Accumulator<PartialLikelihood, double> likelihood;

    // Observation Process
    std::vector<AlleleCounter *> alleleCounters{};
    AlleleCounterAccumulator alleleCountAccumulator;
    ObservationProcessLikelihood<AlleleCounterAccumulator>* observationProcessLikelihood{};

    // Node Transmission Process
    COITransitionProbImpl* coitp{};
    InterTransmissionProbImpl* intp{};
    NodeTransmissionImpl* nodeTransmissionProcess{};

    // Source Transmission Process
    COIProbabilityImpl* coiProb{};
    std::vector<SourceTransmissionImpl *> sourceTransmissionProcessList;

    // Transmission Process
    std::vector<ParentSetImpl *> parentSetList{};
    std::vector<TransmissionProcess *> transmissionProcessList{};

    std::vector<LogTransformer<TransmissionProcess> *> logTransmissionProcessList{};

};


#endif //TRANSMISSION_NETWORKS_APP_MODELONE_H
