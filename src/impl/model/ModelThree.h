////
//// Created by Maxwell Murphy on 6/18/20.
////
//
//#ifndef TRANSMISSION_NETWORKS_APP_MODELTHREE_H
//#define TRANSMISSION_NETWORKS_APP_MODELTHREE_H
//
//#include "core/computation/Accumulator.h"
//
//#include "core/distributions/ZTGeometric.h"
//#include "core/distributions/ZTPoisson.h"
//#include "core/distributions/ZTMultiplicativeBinomial.h"
//
//#include "impl/state/ModelThreeState.h"
//
//#include "model/observation_process/AlleleCounter.h"
//#include "model/observation_process/ObservationProcessLikelihoodv1.h"
//
//#include "model/transmission_process/NetworkBasedTransmissionProcess.h"
//#include "model/transmission_process/node_transmission_process/NoSuperInfectionMutation.h"
//
//#include "model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h"
//
//namespace transmission_nets::impl {
//
//    class ModelThree {
//        using State = ModelThreeState;
//        static constexpr int MAX_COI = 10;
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
//        using InterTransmissionProbImpl = core::distributions::ZTGeometric<MAX_TRANSMISSIONS>;
//        using NodeTransmissionImpl = model::transmission_process::NoSuperInfectionMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>;
//
//        using COIProbabilityImpl = core::distributions::ZTPoisson<MAX_COI>;
//        using SourceTransmissionImpl = model::transmission_process::MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainerImpl, InfectionEvent, MAX_COI>;
//
//        using ParentSetImpl = core::containers::ParentSet<InfectionEvent>;
//        using TransmissionProcess = model::transmission_process::NetworkBasedTransmissionProcess<MAX_PARENTS, NodeTransmissionImpl, SourceTransmissionImpl, InfectionEvent>;
//
////    using BetaPrior = Prior<boost::math::beta_distribution<>, Parameter<float>, int, int>;
////    using GammaPrior = Prior<boost::math::gamma_distribution<>, Parameter<float>, float, float>;
//
//    public:
//        explicit ModelThree(State& state);
//
//        float value();
//
//        bool isDirty();
//
//        State& state;
//        core::computation::Accumulator<core::computation::PartialLikelihood, float> likelihood;
//
//        // Observation Process
//        std::vector<AlleleCounterImpl *> alleleCounters{};
//        AlleleCounterAccumulator alleleCountAccumulator;
//        model::observation_process::ObservationProcessLikelihoodv1<AlleleCounterAccumulator>* observationProcessLikelihood{};
//
//        // Node Transmission Process
//        InterTransmissionProbImpl* intp{};
//        NodeTransmissionImpl* nodeTransmissionProcess{};
//
//        // Source Transmission Process
//        COIProbabilityImpl* coiProb{};
//        std::vector<SourceTransmissionImpl *> sourceTransmissionProcessList;
//
//        // Transmission Process
//        std::vector<TransmissionProcess *> transmissionProcessList{};
//
//    };
//
//
//}
//
//
//#endif //TRANSMISSION_NETWORKS_APP_MODELTHREE_H
