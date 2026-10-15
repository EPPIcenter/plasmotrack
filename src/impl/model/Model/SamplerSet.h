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

        int infection_idx_ = 0;
        for (const auto& infection : state->infections) {

            const bool isSymptomatic = infection->isSymptomatic();
            const std::string infection_id = infection->id();
            const double upperBound = isSymptomatic ? state->symptomaticInfectionDurationDist->value().size() : state->asymptomaticInfectionDurationDist->value().size();

            scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(infection->infectionDuration(), target, 1.0, upperBound, r, 1, .1, 2),
                                        .id = fmt::format("ConstrainedContinuousRandomWalk/infection_duration/{}", infection_id),
                                        .adaptationStart = 20,
                                        .adaptationEnd = 200});
            if (!state->null_model_) {
                scheduler_.registerSampler({.sampler = std::make_unique<genetics::RandomAllelesBitSetSampler3<T, Engine, InfectionEvent, GeneticsImpl, ParentSetImpl, MAX_PARENTS, MAX_COI>>(infection, target->parentSets.at(infection_id), state->latentParents[infection_idx_], target, r),
                                            .id = fmt::format("RandomAllelesBitSetSampler3/latent_genotypes/{}", infection_id)});
                // Move the infection to another source in one step (its latent parent stands for the Source, via its fixed anchor).
                scheduler_.registerSampler({.sampler = std::make_unique<specialized::ParentRelinkSampler<T, Engine, InfectionEvent>>(infection, state->latentParents[infection_idx_], state->allowedRelationships->allowedParents(infection), target, r, 1.0, upperBound, isSymptomatic ? state->symptomaticInfectionDurationDist->value() : state->asymptomaticInfectionDurationDist->value(), 0.01, 0.01, RELINK_TRIES),
                                            .id = fmt::format("ParentRelinkSampler/infection/{}", infection_id)});
                scheduler_.registerSampler({.sampler = std::make_unique<genetics::StrainCountSampler<T, Engine, InfectionEvent>>(infection, target, r, MAX_COI),
                                            .id = fmt::format("StrainCountSampler/latent_genotypes/{}", infection_id)});
                // Move the infection to the other side of a clone in the ordering.
                scheduler_.registerSampler({.sampler = std::make_unique<specialized::CloneOrderSwapSampler<T, Engine, InfectionEvent>>(infection, specialized::CloneOrderSwapSampler<T, Engine, InfectionEvent>::mutualPartners(state->infections, *state->allowedRelationships, infection_idx_), target, r, 1.0, upperBound, isSymptomatic ? state->symptomaticInfectionDurationDist->value() : state->asymptomaticInfectionDurationDist->value()),
                                            .id = fmt::format("CloneOrderSwapSampler/ordering/{}", infection_id)});
                {
                    std::vector<std::shared_ptr<InfectionEvent>> partners;
                    for (const auto& partner : specialized::CloneOrderSwapSampler<T, Engine, InfectionEvent>::mutualPartners(state->infections, *state->allowedRelationships, infection_idx_)) {
                        partners.push_back(partner.infection);
                    }
                    scheduler_.registerSampler({.sampler = std::make_unique<specialized::InformedOrderJumpSampler<T, Engine, InfectionEvent>>(infection, std::move(partners), target, r, 1.0, upperBound, isSymptomatic ? state->symptomaticInfectionDurationDist->value() : state->asymptomaticInfectionDurationDist->value(), INFORMED_JUMP_MAX_TRIES),
                                                .id = fmt::format("InformedOrderJumpSampler/ordering/{}", infection_id)});
                }
                for (const auto& [locus_label, locus] : state->loci) {
                    if (infection->latentGenotype().contains(locus)) {
                        auto latentGenotype = infection->latentGenotype(locus);
                        scheduler_.registerSampler({.sampler = std::make_unique<genetics::RandomAllelesBitSetSampler4<T, Engine, InfectionEvent, GeneticsImpl, LocusImpl>>(infection, state->latentParents[infection_idx_], locus, state->allowedRelationships, target, r, MAX_COI),
                                                    .id = fmt::format("RandomAllelesBitSetSampler4/latent_genotypes/{}/{}", infection_id, locus_label)});
                        // Swaps pay where the observation says little about which allele is present: every Sampling iteration at
                        // loci with no allele observed, every 10th at loci with several, and never at loci with exactly one.
                        const unsigned int observedAlleles = infection->observedGenotype().contains(locus) ? infection->observedGenotype(locus)->value().totalPositiveCount() : 0;
                        if (observedAlleles != 1) {
                            scheduler_.registerSampler({.sampler = std::make_unique<genetics::AlleleSwapSampler<T, Engine, GeneticsImpl>>(infection->latentGenotype(locus), target, r),
                                                    .id = fmt::format("AlleleSwapSampler/latent_genotypes/{}/{}", infection_id, locus_label),
                                                    .updateFrequency = observedAlleles == 0 ? 1 : 10});
                        }
                    }
                }
            }
            infection_idx_++;
        }


        if (!state->null_model_) {

            for (size_t i = 0; i < state->expectedFalseNegatives.size(); ++i) {
                scheduler_.registerSampler({
                        .sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state->expectedFalseNegatives[i], target, 1e-6, .5, r, 1, .1, 2),
                        .id = fmt::format("ConstrainedContinuousRandomWalk/expected_false_negatives/{}", state->infections[i]->id()),
                        .adaptationStart = 20,
                        .adaptationEnd = 200
                });
            }

            for (size_t i = 0; i < state->expectedFalsePositives.size(); ++i) {
                scheduler_.registerSampler({.sampler = std::make_unique<ConstrainedContinuousRandomWalk<T, Engine>>(state->expectedFalsePositives[i], target, 1e-6, .5, r, 1, .1, 2),
                                            .id = fmt::format("ConstrainedContinuousRandomWalk/expected_false_positives/{}", state->infections[i]->id()),
                                            .adaptationStart = 20,
                                            .adaptationEnd = 200});
            }

            for (const auto& [locus_label, locus] : state->loci) {
                // A one-allele locus has a fixed frequency of 1 and nothing to sample.
                if (locus->totalAlleles() < 2) {
                    continue;
                }
                auto alleleFrequencySampler = std::make_unique<SALTSampler<T, Engine>>(state->alleleFrequencies->alleleFrequencies(locus), target, r, 1, .1, 2);
                alleleFrequencySampler->setLowerLimit_(MIN_ALLELE_FREQUENCY);
                scheduler_.registerSampler({.sampler = std::move(alleleFrequencySampler),
                                            .id = fmt::format("SALTSampler/allele_frequencies/{}", locus_label),
                                            .adaptationStart = 20,
                                            .adaptationEnd = 200,
                                            .updateFrequency = SALT_UPDATE_FREQUENCY});
            }
        }

    }
}// namespace transmission_nets::impl::Model

#endif//TRANSMISSION_NETWORKS_APP_SAMPLERSET_H
