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
#include "model/transmission_process/node_transmission_process/LatentParentMarginal.h"
#include "model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h"

#include <filesystem>

namespace transmission_nets::impl::Model {
    // Tries per ParentRelinkSampler update (1: a plain independence move; more: a multiple-try move).
    constexpr unsigned int RELINK_TRIES = 2;
    // Partners tried per InformedOrderJumpSampler update (a random subset when an infection has more).
    constexpr size_t INFORMED_JUMP_MAX_TRIES = 8;
    // Allele-frequency (SALT) updates run every this many Sampling iterations.
    constexpr int SALT_UPDATE_FREQUENCY = 3;
    // Frequencies each SALT update proposes, chosen at random (0: every frequency at the locus, in random order).
    constexpr unsigned int SALT_PROPOSALS_PER_UPDATE = 0;
    // Whether SALT samplers are registered.
    constexpr bool SALT_SAMPLER = true;
    // Proposals per SimplexLogRatioSampler update, which moves all of a locus' frequencies at once (0: not registered).
    constexpr unsigned int FREQUENCY_WALK_PROPOSALS = 0;
    // SimplexLogRatioSampler adapts its scale and covariance until this Sampling iteration.
    constexpr int FREQUENCY_WALK_ADAPTATION_END = 500;
    // Whether JointGeneticsTimeSampler is registered.
    constexpr bool JOINT_GENETICS_TIME_SAMPLER = false;


    static constexpr int MAX_ALLELES       = 128;
    static constexpr int MAX_COI           = 20;
    static constexpr int MAX_PARENTS       = 2;
    static constexpr int MAX_PARENT_SET_SIZE = MAX_PARENTS + 1;
    // static constexpr int MAX_TRANSMISSIONS = 8;
    static constexpr int MAX_STRAINS = 12;
    // Allele frequencies are kept at or above this, at the start and while sampling. A numerical guard only: far below
    // one strain among all sources (about 1 / (2 * infections * mean COI)), so it should never bind.
    static constexpr double MIN_ALLELE_FREQUENCY = 1e-6;

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
    using NodeTransmissionImpl = model::transmission_process::MultinomialTransmissionProcess<MAX_PARENTS, MAX_STRAINS, SourceTransmissionImpl, ParentSetSizeLikelihoodImpl, GeneticsImpl>;

    using OrderingImpl        = core::computation::ObservationTimeDerivedOrdering<InfectionEvent>;
    using ParentSetImpl       = core::computation::OrderDerivedParentSet<InfectionEvent, OrderingImpl>;

    // The latent parent summed out of each infection's transmission term.
    using LatentMarginalImpl = model::transmission_process::LatentParentMarginal<MAX_PARENTS, MAX_COI, MAX_STRAINS, NodeTransmissionImpl, COIProbabilityImpl, AlleleFrequencyContainerImpl, GeneticsImpl, LocusImpl>;

    using TransmissionProcess = model::transmission_process::OrderBasedTransmissionProcessV3<MAX_PARENTS, NodeTransmissionImpl, LatentMarginalImpl, ParentSetSizeLikelihoodImpl, InfectionEvent, ParentSetImpl>;

}// namespace transmission_nets::impl::Model

#endif//TRANSMISSION_NETWORKS_APP_MODEL_CONFIG_H
