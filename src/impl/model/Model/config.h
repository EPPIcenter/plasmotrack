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
    // ParentRelinkSampler weights a candidate parent by exp(scale * (observed alleles it carries - observed alleles it
    // lacks)); 0 averages that score over observed loci instead, which leaves the candidates nearly uniform.
    constexpr double RELINK_CANDIDATE_WEIGHT_SCALE = 1.0;
    // Partners tried per InformedOrderJumpSampler update (a random subset when an infection has more).
    constexpr size_t INFORMED_JUMP_MAX_TRIES = 4;
    // Allele-frequency (SALT) updates run every this many Sampling iterations.
    constexpr int SALT_UPDATE_FREQUENCY = 30;
    // Each allele-frequency sampler, updated every iteration, also makes a Dirichlet move every
    // ALLELE_FREQUENCY_DIRICHLET_EVERY-th update, fitted to its own trajectory over these spans of updates and fixed
    // after the last; its SALT sweep then runs every SALT_UPDATE_FREQUENCY-th update.
    constexpr bool ALLELE_FREQUENCY_DIRICHLET_MOVE = true;
    constexpr int ALLELE_FREQUENCY_DIRICHLET_EVERY = 5;
    inline const std::vector<int> ALLELE_FREQUENCY_FIT_SPANS{5, 45, 100, 200, 400};
    constexpr double ALLELE_FREQUENCY_DIRICHLET_WIDEN = 0.25;
    // A locus of more alleles than this is moved in random blocks of this many; 0 moves the whole locus at once.
    constexpr size_t ALLELE_FREQUENCY_DIRICHLET_BLOCK = 4;
    // CloneOrderSwapSampler also exchanges an infection's infection time with a clone's.
    constexpr bool CLONE_TIME_EXCHANGE = true;
    // Single observed-allele toggles each StrainCountSampler update makes after it adds or removes a strain.
    constexpr unsigned int OBSERVED_ALLELE_TOGGLES = 2;
    // Share of RandomAllelesBitSetSampler3 updates whose parent-set option is uniform; the rest follow the
    // infection's parent-set posterior. 1 makes every option uniform.
    constexpr double PARENT_SET_PROPOSAL_UNIFORM_SHARE = 0.1;
    // Loci RandomAllelesBitSetSampler3 redraws per update, chosen uniformly; 0 redraws every locus.
    constexpr size_t PARENT_SET_PROPOSAL_LOCI = 4;
    // Guided swaps GuidedPairSampler makes per Sampling iteration, each on a pair of infections that can infect each
    // other whose observed genotypes differ at no more than GUIDED_PAIR_MAX_DIFFERENT of the loci observed in both
    // (every such pair when negative), picked uniformly.
    constexpr int GUIDED_PAIR_SWAPS         = 30;
    constexpr int GUIDED_PAIR_MAX_DIFFERENT = 10;
    // With a leak prior, how often the leak takes a step: every this many Sampling iterations. Each step re-evaluates
    // every Transmission term at the proposed leak.
    constexpr int LEAK_UPDATE_FREQUENCY = 1;
    // Hierarchical error model: walks per Sampling iteration on each population mean and shape. They touch only prior
    // terms, so they are cheap.
    constexpr int HIERARCHICAL_RATE_STEPS = 5;
    // Family moves, on the Families of the Mutual pairs: a pair is close when, in the direction where the count is
    // smaller, at most FAMILY_MAX_NOVELTY Loci have an Allele observed in one Infection and not the other (Loci
    // missing in either do not count). Each Sampling iteration makes FAMILY_RELOCATIONS relocations, each moving
    // FAMILY_RELOCATE_MOVERS members to an infection time drawn uniformly over the times their infection duration
    // allows or within FAMILY_RELOCATE_JITTER days of another member's, and FAMILY_PAIR_LOCUS_PROPOSALS pair-locus
    // Gibbs Proposals.
    constexpr int FAMILY_RELOCATIONS             = 5;
    constexpr int FAMILY_PAIR_LOCUS_PROPOSALS    = 20;
    constexpr unsigned int FAMILY_MAX_NOVELTY    = 3;
    constexpr size_t FAMILY_RELOCATE_MOVERS      = 2;
    constexpr double FAMILY_RELOCATE_JITTER      = 10.0;
    // The shortest infection duration, in days, the duration samplers allow.
    constexpr double MIN_INFECTION_DURATION = 1.0;


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
