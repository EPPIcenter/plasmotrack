//
// Created by mmurphy on 10/29/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_MODEL_NINE_CONFIG_H
#define TRANSMISSION_NETWORKS_APP_MODEL_NINE_CONFIG_H

// #include "core/computation/Accumulator.h"
#include "core/computation/ObservationTimeDerivedOrdering.h"
#include "core/computation/PartialLikelihood.h"

#include "core/parameters/Parameter.h"

#include "core/containers/AlleleFrequencyContainer.h"

#include "core/datatypes/Alleles.h"

#include "core/distributions/ZTGeometric.h"
#include "core/distributions/ZTMultiplicativeBinomial.h"
#include "core/distributions/ZTPoisson.h"
#include "core/distributions/pdfs/BetaLogPDF.h"
#include "core/distributions/pdfs/GammaLogPDF.h"
#include "core/distributions/DiscreteDistribution.h"
#include "core/distributions/pdfs/DiscretePDF.h"

#include "core/io/loggers/AbstractLogger.h"
#include "core/io/loggers/FileOutput.h"
#include "core/io/loggers/MultiValueLogger.h"
#include "core/io/loggers/ValueLogger.h"

#include "core/samplers/general/ConstrainedContinuousRandomWalk.h"
#include "core/samplers/general/ConstrainedDiscreteRandomWalk.h"
#include "core/samplers/general/SALTSampler.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler2.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler3.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler4.h"
#include "core/samplers/genetics/SequentialAllelesBitSetSampler.h"
#include "core/samplers/genetics/ZanellaAllelesBitSetSampler.h"
#include "core/samplers/specialized/JointGeneticsTimeSampler.h"
#include "core/samplers/scheduler/RandomizedScheduler.h"

#include "model/observation_process/AlleleCounter.h"
#include "model/observation_process/AlleleCounts.h"
#include "model/observation_process/ObservationProcessLikelihoodv2.h"

#include "model/transmission_process/OrderBasedTransmissionProcessV3.h"
// #include "model/transmission_process/node_transmission_process/SimpleLossMutation.h"
// #include "model/transmission_process/node_transmission_process/SimpleLoss.h"
#include "model/transmission_process/node_transmission_process/MultinomialTransmissionProcess2.h"

#include "model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h"

#include <filesystem>
#include <utility>


namespace transmission_nets::impl::ModelNine {

    static constexpr int MAX_ALLELES       = 64;
    static constexpr int MAX_COI           = 30;
    static constexpr int MAX_PARENTS       = 2;
    static constexpr int MAX_TRANSMISSIONS = 8;
    static constexpr int MAX_STRAINS = 12;

    namespace fs                           = std::filesystem;

    using Likelihood                   = core::computation::Likelihood;
    using LocusImpl                    = core::containers::Locus;
    using GeneticsImpl                 = core::datatypes::AllelesBitSet<MAX_ALLELES>;
    using InfectionEvent               = core::containers::Infection<GeneticsImpl>;
    using AlleleFrequencyImpl          = core::datatypes::Simplex;
    using AlleleFrequencyContainerImpl = core::containers::AlleleFrequencyContainer<AlleleFrequencyImpl>;

    using ObservationProcessImpl   = model::observation_process::ObservationProcessLikelihoodv2<GeneticsImpl>;

    using COIProbabilityImpl     = core::distributions::ZTPoisson<MAX_COI>;
    using SourceTransmissionImpl = model::transmission_process::MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainerImpl, InfectionEvent::GenotypeParameterMap , MAX_COI>;
//    using InterTransmissionProbImpl = core::distributions::ZTGeometric<MAX_TRANSMISSIONS>;
//    using NodeTransmissionImpl      = model::transmission_process::SimpleLoss<MAX_TRANSMISSIONS, MAX_PARENTS, InterTransmissionProbImpl, SourceTransmissionImpl>;
    using NodeTransmissionImpl = model::transmission_process::MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, SourceTransmissionImpl>;


    using OrderingImpl        = core::computation::ObservationTimeDerivedOrdering<InfectionEvent>;
    using ParentSetImpl       = core::computation::OrderDerivedParentSet<InfectionEvent, OrderingImpl>;
    using TransmissionProcess = model::transmission_process::OrderBasedTransmissionProcessV3<MAX_PARENTS, NodeTransmissionImpl, SourceTransmissionImpl, InfectionEvent, ParentSetImpl>;

}// namespace transmission_nets::impl::ModelNine

#endif//TRANSMISSION_NETWORKS_APP_MODEL_NINE_CONFIG_H
