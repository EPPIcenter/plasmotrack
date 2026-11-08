//
// The Sampler set: the samplers attached to one chain, in the fixed order they update in, each with its update
// frequency and adaptation window. Built once per chain from the State and the Model, and the one place that decides
// which samplers a Run uses. Every chain in a Run has the same Sampler set, so sampler ids are the same across chains
// and across Runs from the same input.
//

#ifndef TRANSMISSION_NETWORKS_APP_SAMPLERSET_H
#define TRANSMISSION_NETWORKS_APP_SAMPLERSET_H

#include "config.h"
#include "impl/model/Model/State.h"

#include "core/samplers/AbstractSampler.h"
#include "core/samplers/general/ConstrainedContinuousRandomWalk.h"
#include "core/samplers/general/SALTSampler.h"
#include "core/samplers/genetics/AlleleSwapSampler.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler3.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler4.h"
#include "core/samplers/genetics/StrainCountSampler.h"
#include "core/samplers/scheduler/Scheduler.h"
#include "core/samplers/specialized/CloneOrderSwapSampler.h"
#include "core/samplers/specialized/GuidedPairSampler.h"
#include "core/samplers/specialized/InformedOrderJumpSampler.h"
#include "core/samplers/specialized/ParentRelinkSampler.h"

#include <boost/random.hpp>

#include <fmt/format.h>

#include <memory>
#include <string>
#include <vector>

namespace transmission_nets::impl::Model {

    template<typename T, typename Engine = boost::random::mt19937>
    class SamplerSet {
    public:
        SamplerSet(std::shared_ptr<State> state, std::shared_ptr<T> target, std::shared_ptr<Engine> r);

        // One Sampling iteration: every sampler due this iteration updates, and those inside their adaptation window
        // also adapt, unless `adapt` is false.
        void step(const bool adapt = true) {
            // The allele-frequency samplers' Dirichlet fits follow their own trajectories; with fewer updates than the
            // fitting spans in burn-in they would still be refitting while sampling. Fix them when adaptation stops.
            if (not adapt and not fitsFrozen_) {
                for (auto* salt : saltSamplers_) {
                    salt->freezeFit();
                }
                fitsFrozen_ = true;
            }
            scheduler_.step(adapt);
        }

        // Sampler ids, `<type>/<target>`, in update order.
        [[nodiscard]] const std::vector<std::string>& ids() const noexcept {
            return scheduler_.samplerIds();
        }

        // Throws std::out_of_range for an id that isn't registered.
        [[nodiscard]] core::samplers::AbstractSampler& sampler(const std::string& id) const {
            return scheduler_.sampler(id);
        }

        // What this set carries between Sampling iterations, by sampler in ids() order; see SamplerSetRecord.
        [[nodiscard]] core::samplers::SamplerSetRecord record() const {
            return scheduler_.record();
        }

        // Continues from `saved`, as when resuming from a Checkpoint. Throws std::invalid_argument for a record of
        // another Sampler set.
        void restore(const core::samplers::SamplerSetRecord& saved) {
            scheduler_.restore(saved);
        }

    private:
        core::samplers::Scheduler scheduler_;
        std::vector<core::samplers::SALTSampler<T, Engine>*> saltSamplers_{};
        bool fitsFrozen_ = false;
    };

    template<typename T, typename Engine>
    SamplerSet<T, Engine>::SamplerSet(std::shared_ptr<State> state, std::shared_ptr<T> target, std::shared_ptr<Engine> r) {
        using namespace core::samplers;

        scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state->meanCOI, target, 1.0, 20, r, 1, .1, 1),
                                    .id = "ConstrainedContinuousRandomWalk/mean_coi",
                                    .adaptationStart = 20,
                                    .adaptationEnd = 200});

        scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state->meanStrainsTransmitted, target, 1.0, 20, r, 1, .1, 1),
                                    .id = "ConstrainedContinuousRandomWalk/mean_strains_tx",
                                    .adaptationStart = 20,
                                    .adaptationEnd = 200});

        scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state->parentSetSizeProb, target, 0, 1, r, 1, .1, 1),
                                    .id = "ConstrainedContinuousRandomWalk/parent_set_size_prob",
                                    .adaptationStart = 20,
                                    .adaptationEnd = 200});

        // An Infection's longest duration: the length of its duration prior.
        const auto maxDurationOf = [&](const InfectionEvent& infection) {
            return static_cast<double>((infection.isSymptomatic() ? state->symptomaticInfectionDurationDist : state->asymptomaticInfectionDurationDist)->value().size());
        };

        int infectionIndex = 0;
        for (const auto& infection : state->infections) {

            const bool isSymptomatic = infection->isSymptomatic();
            const std::string infection_id = infection->id();
            const auto durations = isSymptomatic ? state->symptomaticInfectionDurationDist->value() : state->asymptomaticInfectionDurationDist->value();
            const double upperBound = durations.size();

            scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(infection->infectionDuration(), target, 1.0, upperBound, r, 1, .1, 2),
                                        .id = fmt::format("ConstrainedContinuousRandomWalk/infection_duration/{}", infection_id),
                                        .adaptationStart = 20,
                                        .adaptationEnd = 200});
            if (!state->nullModel) {
                genetics::ParentSetLliks parentSetLliks = nullptr;
                if constexpr (requires { target->transmissionProcessList; }) {
                    parentSetLliks = [term = target->transmissionProcessList[infectionIndex]](std::vector<double>& out) { return term->parentSetLliks(out); };
                }
                scheduler_.registerSampler({.sampler = std::make_unique<genetics::RandomAllelesBitSetSampler3<T, Engine, InfectionEvent, OrderingImpl, MAX_PARENTS, MAX_COI>>(infection, state->infectionEventOrdering, state->latentParents[infectionIndex], target, r, 0.01, 0.01, state->expectedFalsePositives[infectionIndex], state->expectedFalseNegatives[infectionIndex], std::move(parentSetLliks), PARENT_SET_PROPOSAL_UNIFORM_SHARE, PARENT_SET_PROPOSAL_LOCI),
                                            .id = fmt::format("RandomAllelesBitSetSampler3/latent_genotypes/{}", infection_id)});
                // Move the infection to another of its candidate parents in one step (its latent parent stands for the Source, via its fixed anchor).
                // The 0.01 rates apply only if the infection's sampled error counts are absent; those counts replace them.
                scheduler_.registerSampler({.sampler = std::make_unique<specialized::ParentRelinkSampler<T, Engine, InfectionEvent>>(infection, state->latentParents[infectionIndex], state->allowedRelationships->permittedParents(infection), target, r, 1.0, upperBound, durations, 0.01, 0.01, RELINK_TRIES, state->expectedFalsePositives[infectionIndex], state->expectedFalseNegatives[infectionIndex], RELINK_CANDIDATE_WEIGHT_SCALE),
                                            .id = fmt::format("ParentRelinkSampler/infection/{}", infection_id)});
                auto strainCountSampler = std::make_unique<genetics::StrainCountSampler<T, Engine, InfectionEvent>>(infection, target, r, MAX_COI);
                strainCountSampler->enableObservedToggles(OBSERVED_ALLELE_TOGGLES);
                scheduler_.registerSampler({.sampler = std::move(strainCountSampler),
                                            .id = fmt::format("StrainCountSampler/latent_genotypes/{}", infection_id)});
                // Move the infection to the other side of a clone in the ordering.
                using CloneOrderSwap = specialized::CloneOrderSwapSampler<T, Engine, InfectionEvent>;
                auto mutualPartners  = CloneOrderSwap::mutualPartners(state->infections, *state->allowedRelationships, infectionIndex);
                std::vector<std::shared_ptr<InfectionEvent>> partners;
                for (auto& partner : mutualPartners) {
                    partners.push_back(partner.infection);
                    partner.maxDuration = maxDurationOf(*partner.infection);
                }
                auto cloneOrderSwap = std::make_unique<CloneOrderSwap>(infection, std::move(mutualPartners), target, r, 1.0, upperBound, durations);
                if constexpr (CLONE_TIME_EXCHANGE) {
                    cloneOrderSwap->enableTimeExchange();
                }
                scheduler_.registerSampler({.sampler = std::move(cloneOrderSwap),
                                            .id = fmt::format("CloneOrderSwapSampler/ordering/{}", infection_id)});
                scheduler_.registerSampler({.sampler = std::make_unique<specialized::InformedOrderJumpSampler<T, Engine, InfectionEvent>>(infection, std::move(partners), target, r, 1.0, upperBound, durations, INFORMED_JUMP_MAX_TRIES),
                                            .id = fmt::format("InformedOrderJumpSampler/ordering/{}", infection_id)});
                for (const auto& [locus_label, locus] : state->loci) {
                    if (infection->latentGenotype().contains(locus)) {
                        // Same schedule as the allele swap below. A locus with one observed allele still has a sampler, every
                        // 10th iteration: dropping it would change the Sampler set and break a Checkpoint. The path walk
                        // reads this infection's sampled error counts.
                        const unsigned int observedAlleles = infection->observedGenotype().contains(locus) ? infection->observedGenotype(locus)->value().totalPositiveCount() : 0;
                        scheduler_.registerSampler({.sampler = std::make_unique<genetics::RandomAllelesBitSetSampler4<T, Engine, InfectionEvent, LocusImpl>>(infection, state->latentParents[infectionIndex], locus, state->allowedRelationships, target, r, MAX_COI, state->expectedFalsePositives[infectionIndex], state->expectedFalseNegatives[infectionIndex]),
                                                    .id = fmt::format("RandomAllelesBitSetSampler4/latent_genotypes/{}/{}", infection_id, locus_label),
                                                    .updateFrequency = observedAlleles == 0 ? 1 : 10});
                        // Swaps pay where the observation says little about which allele is present: every Sampling iteration at
                        // loci with no allele observed, every 10th at loci with several, and never at loci with exactly one.
                        if (observedAlleles != 1) {
                            scheduler_.registerSampler({.sampler = std::make_unique<genetics::AlleleSwapSampler<T, Engine>>(infection->latentGenotype(locus), target, r),
                                                    .id = fmt::format("AlleleSwapSampler/latent_genotypes/{}/{}", infection_id, locus_label),
                                                    .updateFrequency = observedAlleles == 0 ? 1 : 10});
                        }
                    }
                }
            }
            infectionIndex++;
        }


        if (!state->nullModel) {

            // One walk per Infection on each observation rate, registered rate by rate.
            const auto walkEachRate = [&](const auto& rates, const std::string& name) {
                for (size_t i = 0; i < rates.size(); ++i) {
                    scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(rates[i], target, 1e-6, .5, r, 1, .1, 2),
                                                .id = fmt::format("ConstrainedContinuousRandomWalk/{}/{}", name, state->infections[i]->id()),
                                                .adaptationStart = 20,
                                                .adaptationEnd = 200});
                }
            };
            walkEachRate(state->expectedFalseNegatives, "expected_false_negatives");
            walkEachRate(state->expectedFalsePositives, "expected_false_positives");

            for (const auto& [locus_label, locus] : state->loci) {
                // A one-allele locus has a fixed frequency of 1 and nothing to sample.
                if (locus->totalAlleles() < 2) {
                    continue;
                }
                auto alleleFrequencySampler = std::make_unique<SALTSampler<T, Engine>>(state->alleleFrequencies->alleleFrequencies(locus), target, r, 1, .1, 2);
                alleleFrequencySampler->setLowerLimit_(MIN_ALLELE_FREQUENCY);
                if constexpr (ALLELE_FREQUENCY_DIRICHLET_MOVE) {
                    alleleFrequencySampler->enableDirichletMove(ALLELE_FREQUENCY_FIT_SPANS, ALLELE_FREQUENCY_DIRICHLET_WIDEN, SALT_UPDATE_FREQUENCY, ALLELE_FREQUENCY_DIRICHLET_EVERY,
                                                                ALLELE_FREQUENCY_DIRICHLET_BLOCK);
                }
                saltSamplers_.push_back(alleleFrequencySampler.get());
                scheduler_.registerSampler({.sampler = std::move(alleleFrequencySampler),
                                            .id = fmt::format("SALTSampler/allele_frequencies/{}", locus_label),
                                            .adaptationStart = 20,
                                            .adaptationEnd = 200,
                                            .updateFrequency = ALLELE_FREQUENCY_DIRICHLET_MOVE ? 1 : SALT_UPDATE_FREQUENCY});
            }

            // Last, so that the samplers before it keep their positions in a Checkpoint's record.
            using GuidedPairs = specialized::GuidedPairSampler<T, Engine, InfectionEvent, LocusImpl>;
            auto maxDuration  = [&](const size_t i) { return maxDurationOf(*state->infections[i]); };
            scheduler_.registerSampler({.sampler = std::make_unique<GuidedPairs>(GuidedPairs::mutualPairs(state->infections, *state->allowedRelationships, maxDuration), target, r, 1.0, MAX_COI,
                                                                                 typename GuidedPairs::Settings{.moves = GUIDED_PAIR_SWAPS, .maxDifferent = GUIDED_PAIR_MAX_DIFFERENT}),
                                        .id      = "GuidedPairSampler/pairs"});
        }

    }
}// namespace transmission_nets::impl::Model

#endif//TRANSMISSION_NETWORKS_APP_SAMPLERSET_H
