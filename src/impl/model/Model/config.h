//
// Created by mmurphy on 10/29/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_MODEL_CONFIG_H
#define TRANSMISSION_NETWORKS_APP_MODEL_CONFIG_H

#include "core/computation/ObservationTimeDerivedOrdering.h"
#include "core/computation/PartialLikelihood.h"

#include "core/containers/AlleleFrequencyContainer.h"

#include "core/datatypes/Alleles.h"


#include "core/distributions/ZTPoisson.h"
#include "core/distributions/ZTGeometric.h"

#include "model/observation_process/ObservationProcessLikelihoodv2.h"

#include "model/transmission_process/OrderBasedTransmissionProcessV3.h"
#include "model/transmission_process/node_transmission_process/MultinomialTransmissionProcess.h"
#include "model/transmission_process/node_transmission_process/LatentParentMarginal.h"
#include <filesystem>

namespace transmission_nets::impl::Model {
    // Tries per ParentRelinkSampler update (1: a plain independence move; more: a multiple-try move).
    constexpr unsigned int RELINK_TRIES = 2;
    // Partners tried per InformedOrderJumpSampler update (a random subset when an infection has more).
    constexpr size_t INFORMED_JUMP_MAX_TRIES = 8;
    // Allele-frequency (SALT) updates run every this many Sampling iterations.
    constexpr int SALT_UPDATE_FREQUENCY = 3;


    static constexpr int MAX_ALLELES       = 128;
    static constexpr int MAX_COI           = 20;
    static constexpr int MAX_PARENTS       = 2;
    static constexpr int MAX_PARENT_SET_SIZE = MAX_PARENTS + 1;
    static constexpr int MAX_STRAINS = 12;
    // Allele frequencies are kept at or above this, at the start and while sampling. A numerical guard only: far below
    // one strain among all sources (about 1 / (2 * infections * mean COI)), so it should never bind.
    static constexpr double MIN_ALLELE_FREQUENCY = 1e-6;

    namespace fs                           = std::filesystem;

    using Likelihood                   = core::computation::Likelihood;
    using LocusImpl                    = core::containers::Locus;
    using GeneticsImpl                 = core::datatypes::AllelesBitSet<MAX_ALLELES>;
    using InfectionEvent               = core::containers::Infection<GeneticsImpl>;
    using AlleleFrequencyImpl          = core::datatypes::Simplex;
    using AlleleFrequencyContainerImpl = core::containers::AlleleFrequencyContainer<AlleleFrequencyImpl>;

    using ObservationProcessImpl   = model::observation_process::ObservationProcessLikelihoodv2<GeneticsImpl>;

    using COIProbabilityImpl     = core::distributions::ZTPoisson<MAX_COI>;
    using ParentSetSizeLikelihoodImpl = core::distributions::ZTGeometric<MAX_PARENT_SET_SIZE>;
    using NodeTransmissionImpl = model::transmission_process::MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, ParentSetSizeLikelihoodImpl, GeneticsImpl>;

    using OrderingImpl        = core::computation::ObservationTimeDerivedOrdering<InfectionEvent>;
    using ParentSetImpl       = core::computation::OrderDerivedParentSet<InfectionEvent, OrderingImpl>;

    // The latent parent summed out of each infection's transmission term.
    using LatentMarginalImpl = model::transmission_process::LatentParentMarginal<MAX_PARENTS, MAX_COI, MAX_STRAINS, NodeTransmissionImpl, COIProbabilityImpl, AlleleFrequencyContainerImpl, GeneticsImpl, LocusImpl>;

    using TransmissionProcess = model::transmission_process::OrderBasedTransmissionProcessV3<MAX_PARENTS, NodeTransmissionImpl, LatentMarginalImpl, ParentSetSizeLikelihoodImpl, InfectionEvent, ParentSetImpl>;

}// namespace transmission_nets::impl::Model

#endif//TRANSMISSION_NETWORKS_APP_MODEL_CONFIG_H
