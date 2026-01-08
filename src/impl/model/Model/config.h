//
// Created by mmurphy on 10/29/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_MODEL_CONFIG_H
#define TRANSMISSION_NETWORKS_APP_MODEL_CONFIG_H

#include "core/computation/ObservationTimeDerivedOrdering.h"
#include "core/computation/PartialLikelihood.h"

#include "core/containers/AlleleFrequencyContainer.h"

#include "core/datatypes/Alleles.h"
// #include "core/datatypes/SparseAlleleSet.h"


#include "core/distributions/ZTPoisson.h"
#include "core/distributions/ZTGeometric.h"

#include "model/observation_process/ObservationProcessLikelihoodv2.h"

#include "model/transmission_process/OrderBasedTransmissionProcessV3.h"
#include "model/transmission_process/node_transmission_process/MultinomialTransmissionProcess.h"
#include "model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h"

#include <filesystem>

namespace transmission_nets::impl::Model {

    static constexpr int MAX_ALLELES       = 128;
    static constexpr int MAX_COI           = 20;
    static constexpr int MAX_PARENTS       = 2;
    static constexpr int MAX_PARENT_SET_SIZE = MAX_PARENTS + 1;
    // static constexpr int MAX_TRANSMISSIONS = 8;
    static constexpr int MAX_STRAINS = 12;

    namespace fs                           = std::filesystem;

    using Likelihood                   = core::computation::Likelihood;
    using LocusImpl                    = core::containers::Locus;
    using GeneticsImpl                 = core::datatypes::AllelesBitSet<MAX_ALLELES>;
    // using GeneticsImpl                 = core::datatypes::SparseAlleleSet;
    using InfectionEvent               = core::containers::Infection<GeneticsImpl>;
    using AlleleFrequencyImpl          = core::datatypes::Simplex;
    using AlleleFrequencyContainerImpl = core::containers::AlleleFrequencyContainer<AlleleFrequencyImpl>;

    using ObservationProcessImpl   = model::observation_process::ObservationProcessLikelihoodv2<GeneticsImpl>;

    using COIProbabilityImpl     = core::distributions::ZTPoisson<MAX_COI>;
    using ParentSetSizeLikelihoodImpl = core::distributions::ZTGeometric<MAX_PARENT_SET_SIZE>;
    using SourceTransmissionImpl = model::transmission_process::MultinomialSourceTransmissionProcess<COIProbabilityImpl, AlleleFrequencyContainerImpl, InfectionEvent::GenotypeParameterMap , MAX_COI>;
//    using InterTransmissionProbImpl = core::distributions::ZTGeometric<MAX_TRANSMISSIONS>;
//    using NodeTransmissionImpl      = model::transmission_process::SimpleLoss<MAX_TRANSMISSIONS, MAX_PARENTS, InterTransmissionProbImpl, SourceTransmissionImpl>;
    using NodeTransmissionImpl = model::transmission_process::MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, SourceTransmissionImpl, ParentSetSizeLikelihoodImpl>;

    using OrderingImpl        = core::computation::ObservationTimeDerivedOrdering<InfectionEvent>;
    using ParentSetImpl       = core::computation::OrderDerivedParentSet<InfectionEvent, OrderingImpl>;

    using TransmissionProcess = model::transmission_process::OrderBasedTransmissionProcessV3<MAX_PARENTS, NodeTransmissionImpl, SourceTransmissionImpl, ParentSetSizeLikelihoodImpl, InfectionEvent, ParentSetImpl>;

}// namespace transmission_nets::impl::Model

#endif//TRANSMISSION_NETWORKS_APP_MODEL_CONFIG_H
