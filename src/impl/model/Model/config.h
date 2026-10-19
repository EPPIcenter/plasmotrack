//
// Created by mmurphy on 10/29/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_MODEL_CONFIG_H
#define TRANSMISSION_NETWORKS_APP_MODEL_CONFIG_H

#include "core/computation/ObservationTimeDerivedOrdering.h"
#include "core/computation/PartialLikelihood.h"

#include "core/containers/AlleleFrequencyContainer.h"

#include "core/datatypes/Genotype.h"

#include <vector>


#include "core/distributions/ZTPoisson.h"
#include "core/distributions/ZTGeometric.h"

#include "model/observation_process/ObservationProcessLikelihoodv2.h"

#include "model/transmission_process/TransmissionTerm.h"
#include "model/transmission_process/TransmissionMemo.h"
#include <filesystem>

namespace transmission_nets::impl::Model {
    // Tries per ParentRelinkSampler update (1: a plain independence move; more: a multiple-try move).
    constexpr unsigned int RELINK_TRIES = 1;
    // ParentRelinkSampler weights a source by exp(scale * (observed alleles it carries - observed alleles it lacks));
    // 0 averages that score over observed loci instead, which leaves the sources nearly uniform.
    constexpr double RELINK_SOURCE_WEIGHT_SCALE = 1.0;
    // Partners tried per InformedOrderJumpSampler update (a random subset when an infection has more).
    constexpr size_t INFORMED_JUMP_MAX_TRIES = 4;
    // Allele-frequency (SALT) updates run every this many Sampling iterations.
    constexpr int SALT_UPDATE_FREQUENCY = 30;
    // Each allele-frequency sampler, updated every iteration, also makes a Dirichlet move every
    // ALLELE_FREQUENCY_DIRICHLET_EVERY-th update, fitted to its own trajectory over these windows of updates and fixed
    // after the last; its SALT sweep then runs every SALT_UPDATE_FREQUENCY-th update.
    constexpr bool ALLELE_FREQUENCY_DIRICHLET_MOVE = true;
    constexpr int ALLELE_FREQUENCY_DIRICHLET_EVERY = 2;
    inline const std::vector<int> ALLELE_FREQUENCY_FIT_WINDOWS{20, 50, 100, 200, 400};
    constexpr double ALLELE_FREQUENCY_DIRICHLET_WIDEN = 0.25;
    // A locus of more alleles than this is moved in random blocks of this many; 0 moves the whole locus at once.
    constexpr size_t ALLELE_FREQUENCY_DIRICHLET_BLOCK = 4;
    // CloneOrderSwapSampler also exchanges an infection's infection time with a clone's.
    constexpr bool CLONE_TIME_EXCHANGE = true;
    // Single observed-allele toggles each StrainCountSampler update makes after its strain-count move.
    constexpr unsigned int OBSERVED_ALLELE_TOGGLES = 8;
    // Share of RandomAllelesBitSetSampler3 updates whose parent-set option is uniform; the rest follow the
    // infection's parent-set posterior. 1 makes every option uniform.
    constexpr double PARENT_SET_PROPOSAL_UNIFORM_SHARE = 0.1;
    // Loci RandomAllelesBitSetSampler3 redraws per update, chosen uniformly; 0 redraws every locus.
    constexpr size_t PARENT_SET_PROPOSAL_LOCI = 4;
    // Guided swaps GuidedPairSampler makes per Sampling iteration, each on a pair of infections that can infect each
    // other whose observed genotypes differ at no more than GUIDED_PAIR_MAX_DIFFERENT of the loci observed in both
    // (every such pair when negative).
    constexpr int GUIDED_PAIR_SWAPS         = 30;
    constexpr int GUIDED_PAIR_MAX_DIFFERENT = 10;


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
    using core::datatypes::Genotype;
    using InfectionEvent               = core::containers::Infection<>;
    using AlleleFrequencyImpl          = core::datatypes::Simplex;
    using AlleleFrequencyContainerImpl = core::containers::AlleleFrequencyContainer<AlleleFrequencyImpl>;

    using ObservationProcessImpl   = model::observation_process::ObservationProcessLikelihoodv2;

    using COIProbabilityImpl     = core::distributions::ZTPoisson<MAX_COI>;
    using ParentSetSizeLikelihoodImpl = core::distributions::ZTGeometric<MAX_PARENT_SET_SIZE>;
    using OrderingImpl        = core::computation::ObservationTimeDerivedOrdering<InfectionEvent>;

    // The per-chain memory of Transmission kernel results, with the latent parent summed out of each term through it.
    using TransmissionMemoImpl = model::transmission_process::TransmissionMemo<MAX_PARENTS, MAX_STRAINS, MAX_COI, InfectionEvent, LocusImpl, AlleleFrequencyContainerImpl>;

    using TransmissionProcess = model::transmission_process::TransmissionTerm<MAX_PARENTS, TransmissionMemoImpl, COIProbabilityImpl, ParentSetSizeLikelihoodImpl, InfectionEvent, OrderingImpl>;

}// namespace transmission_nets::impl::Model

#endif//TRANSMISSION_NETWORKS_APP_MODEL_CONFIG_H
