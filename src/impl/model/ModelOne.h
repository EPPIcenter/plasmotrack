////
//// Created by Maxwell Murphy on 4/19/20.
////
//
//#ifndef TRANSMISSION_NETWORKS_APP_MODELONE_H
//#define TRANSMISSION_NETWORKS_APP_MODELONE_H
//
//#include <utility>
//#include <boost/math/distributions/beta.hpp>
//#include <boost/math/distributions/gamma.hpp>
//
//#include "core/computation/PartialLikelihood.h"
//#include "core/computation/Accumulator.h"
//#include "core/computation/transformers/LogTransformer.h"
//
//#include "core/datatypes/Alleles.h"
//
//#include "core/distributions/ZTGeometric.h"
//#include "core/distributions/ZTMultiplicativeBinomial.h"
//
//#include "core/priors/Prior.h"
//
//#include "impl/state/ModelOneState.h"
//
//#include "model/observation_process/AlleleCounter.h"
//#include "model/observation_process/AlleleCounts.h"
//#include "model/observation_process/ObservationProcessLikelihoodv1.h"
//
//#include "model/transmission_process/OrderBasedTransmissionProcess.h"
//#include "model/transmission_process/node_transmission_process/NoSuperInfectionNoMutation.h"
//
//#include "model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h"
//
//namespace transmission_nets::impl {
//
//    using Likelihood = core::computation::Likelihood;
//    class ModelOne {
//        using State = ModelOneState;
//        static constexpr int MAX_COI = 32;
//        static constexpr int MAX_PARENTS = 1;
//        static constexpr int MAX_TRANSMISSIONS = 5;
//
//        using GeneticsImpl = State::GeneticsImpl;
//        using InfectionEvent = State::InfectionEvent;
//        using AlleleFrequencyContainerImpl = State::AlleleFrequencyContainerImpl;
//
//        using AlleleCounterImpl = model::observation_process::AlleleCounter<GeneticsImpl>;
//        using AlleleCounterAccumulator = core::computation::Accumulator<AlleleCounterImpl, model::observation_process::AlleleCounts>;
//
//        using COITransitionProbImpl = core::distributions::ZTMultiplicativeBinomial<MAX_COI>;
//        using InterTransmissionProbImpl = core::distributions::ZTGeometric<MAX_TRANSMISSIONS>;
//        using NodeTransmissionImpl = model::transmission_process::NoSuperInfectionNoMutation<MAX_COI, MAX_TRANSMISSIONS, COITransitionProbImpl, InterTransmissionProbImpl>;
//
//        using COIProbabilityImpl = core::distributions::ZTGeometric<MAX_COI>;
//        using SourceTransmissionImpl = model::transmission_process::MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainerImpl, InfectionEvent, MAX_COI>;
//
//        using OrderingImpl = core::parameters::Ordering<InfectionEvent>;
//        using ParentSetImpl = core::computation::OrderDerivedParentSet<InfectionEvent, OrderingImpl>;
//        using TransmissionProcess = model::transmission_process::OrderBasedTransmissionProcess<MAX_PARENTS, NodeTransmissionImpl, SourceTransmissionImpl, InfectionEvent, ParentSetImpl>;
//
//        using BetaPrior = core::priors::Prior<boost::math::beta_distribution<>, core::parameters::Parameter<double>, int, int>;
//        using GammaPrior = core::priors::Prior<boost::math::gamma_distribution<>, core::parameters::Parameter<double>, double, double>;
//
//    public:
//        explicit ModelOne(State& state);
//
//        Likelihood value();
//
//        bool isDirty();
//
//        State& state;
//        core::computation::Accumulator<core::computation::PartialLikelihood, Likelihood> likelihood;
//
//        // Observation Process
//        std::vector<AlleleCounterImpl *> alleleCounters{};
//        AlleleCounterAccumulator alleleCountAccumulator;
//        model::observation_process::ObservationProcessLikelihoodv1<AlleleCounterAccumulator>* observationProcessLikelihood{};
//
//        // Node Transmission Process
//        COITransitionProbImpl* coitp{};
//        InterTransmissionProbImpl* intp{};
//        NodeTransmissionImpl* nodeTransmissionProcess{};
//
//        // Source Transmission Process
//        COIProbabilityImpl* coiProb{};
//        std::vector<SourceTransmissionImpl *> sourceTransmissionProcessList;
//
//        // Transmission Process
//        std::vector<ParentSetImpl *> parentSetList{};
//        std::vector<TransmissionProcess *> transmissionProcessList{};
//
//        std::vector<core::computation::LogTransformer<TransmissionProcess> *> logTransmissionProcessList{};
//
//    };
//
//}
//
//
//
//
//#endif //TRANSMISSION_NETWORKS_APP_MODELONE_H
