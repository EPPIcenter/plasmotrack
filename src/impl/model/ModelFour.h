//
// Created by Maxwell Murphy on 12/10/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MODELFOUR_H
#define TRANSMISSION_NETWORKS_APP_MODELFOUR_H

#include <utility>
#include <filesystem>

#include "core/computation/Accumulator.h"
#include "core/computation/PartialLikelihood.h"
#include "core/computation/Transformers/LogTransformer.h"

#include "core/containers/AlleleFrequencyContainer.h"

#include "core/datatypes/Alleles.h"

#include "core/distributions/ZTGeometric.h"
#include "core/distributions/ZTPoisson.h"
#include "core/distributions/ZTMultiplicativeBinomial.h"

#include "core/io/loggers/AbstractLogger.h"

#include "model/observation_process/AlleleCounter.h"
#include "model/observation_process/AlleleCounts.h"
#include "model/observation_process/ObservationProcessLikelihood.h"

#include "model/transmission_process/OrderBasedTransmissionProcess.h"
#include "model/transmission_process/node_transmission_process/NoSuperInfectionMutation.h"

#include "model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h"


namespace transmission_nets::impl::ModelFour {
    static constexpr int MAX_ALLELES = 32;
    static constexpr int MAX_COI = 8;
    static constexpr int MAX_PARENTS = 1;
    static constexpr int MAX_TRANSMISSIONS = 5;
    namespace fs = std::filesystem;

    using Likelihood = core::computation::Likelihood;
    using LocusImpl = core::containers::Locus;
    using GeneticsImpl = core::datatypes::AllelesBitSet<MAX_ALLELES>;
    using InfectionEvent = core::containers::Infection<GeneticsImpl, LocusImpl>;
    using AlleleFrequencyImpl = core::datatypes::Simplex;
    using AlleleFrequencyContainerImpl = core::containers::AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>;

    using AlleleCounterImpl = model::observation_process::AlleleCounter<GeneticsImpl>;
    using AlleleCounterAccumulator = core::computation::Accumulator<AlleleCounterImpl, model::observation_process::AlleleCounts>;

    using InterTransmissionProbImpl = core::distributions::ZTGeometric<MAX_TRANSMISSIONS>;
    using NodeTransmissionImpl = model::transmission_process::NoSuperInfectionMutation<MAX_TRANSMISSIONS, InterTransmissionProbImpl>;

    using COIProbabilityImpl = core::distributions::ZTPoisson<MAX_COI>;
    using SourceTransmissionImpl = model::transmission_process::MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainerImpl, InfectionEvent, MAX_COI>;

    using ParentSetImpl = core::computation::OrderDerivedParentSet<InfectionEvent>;
    using TransmissionProcess = model::transmission_process::OrderBasedTransmissionProcess<MAX_PARENTS, NodeTransmissionImpl, SourceTransmissionImpl, InfectionEvent>;

    struct State {
        State(std::map<std::string, LocusImpl *> &loci, std::vector<InfectionEvent *> &infections, std::map<InfectionEvent *, std::vector<InfectionEvent *>> &disallowedParents);
        std::map<std::string, LocusImpl *> loci{};
        std::vector<InfectionEvent *> infections{};
        std::map<InfectionEvent *, std::vector<InfectionEvent *>> disallowedParents{};

        AlleleFrequencyContainerImpl alleleFrequencies;

        // Network Structure
        core::parameters::Ordering<InfectionEvent> infectionEventOrdering;

        // Observation Process
        std::vector<core::parameters::Parameter<double>> observationFalsePositiveRates{};
        std::vector<core::parameters::Parameter<double>> observationFalseNegativeRates{};

        // Node Transmission Process
        core::parameters::Parameter<double> geometricGenerationProb;
        core::parameters::Parameter<double> lossProb;
        core::parameters::Parameter<double> mutationProb;

        // Source Transmission Process
        core::parameters::Parameter<double> meanCOI;
    };

    struct Model {
        Model(std::map<std::string, LocusImpl *> &loci, std::vector<InfectionEvent *> &infections, std::map<InfectionEvent *, std::vector<InfectionEvent *>> &disallowedParents);
        Likelihood value();
        bool isDirty();

        State state;

        core::computation::Accumulator<core::computation::PartialLikelihood, Likelihood> likelihood;

        // Observation Process
        std::vector<AlleleCounterImpl *> alleleCounters{};
        std::vector<AlleleCounterAccumulator *> alleleCountAccumulators{};
        model::observation_process::ObservationProcessLikelihood<AlleleCounterAccumulator> *observationProcessLikelihood{};

        // Node Transmission Process
        InterTransmissionProbImpl *intp{};
        NodeTransmissionImpl *nodeTransmissionProcess{};

        // Source Transmission Process
        COIProbabilityImpl *coiProb{};
        std::vector<SourceTransmissionImpl *> sourceTransmissionProcessList;

        // Transmission Process
        std::vector<ParentSetImpl *> parentSetList{};
        std::vector<TransmissionProcess *> transmissionProcessList{};
    };

    struct ParameterLogger {
        ParameterLogger(Model& model, fs::path rootPath);

        void logParameters() const;

        Model& model;
        fs::path rootPath;
        fs::path paramOutputFolder;
        fs::path statOutputFolder;
        std::vector<core::io::AbstractLogger*> loggers{};

    };
}

#endif//TRANSMISSION_NETWORKS_APP_MODELFOUR_H
