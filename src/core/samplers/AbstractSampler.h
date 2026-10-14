//
// A sampler: one Metropolis-Hastings move on some of a chain's parameters, run once per Sampling iteration it is due.
// Each update is one Proposal (or none, when the move has nothing to propose), which counts into the sampler's counts.
//

#ifndef TRANSMISSION_NETWORKS_APP_ABSTRACTSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_ABSTRACTSAMPLER_H

#include "core/computation/PartialLikelihood.h"
#include "core/samplers/Proposal.h"
#include "core/samplers/SamplerCounts.h"

#include <stdexcept>
#include <string>

namespace transmission_nets::core::samplers {
    using Likelihood = core::computation::Likelihood;

    class AbstractSampler {
    public:
        virtual ~AbstractSampler()     = default;
        virtual void update() noexcept = 0;
        virtual void adapt() noexcept {};

        [[nodiscard]] SamplerCounts counts() const noexcept { return counts_; }

        // Continue counting from `counts`, as when resuming from a checkpoint.
        void restoreCounts(const SamplerCounts counts) noexcept { counts_ = counts; }

        [[nodiscard]] virtual AdaptationState adaptationState() const { return {}; }

        // Continue adapting from `state`, as when resuming from a checkpoint. Throws std::invalid_argument when
        // `state` has a different number of scales than this sampler.
        virtual void restoreAdaptationState(const AdaptationState& state) {
            requireScales(state, 0);
        }

    protected:
        static void requireScales(const AdaptationState& state, const size_t expected) {
            if (state.scales.size() != expected) {
                throw std::invalid_argument("adaptation state has " + std::to_string(state.scales.size()) + " proposal scales, sampler has " + std::to_string(expected));
            }
        }

        // Opens a Proposal on `target`, counting into this sampler's counts when it settles.
        template<typename Target, typename Engine>
        [[nodiscard]] Proposal<Target, Engine> propose(Target& target, Engine& rng) {
            return Proposal<Target, Engine>(target, rng, counts_);
        }

    private:
        SamplerCounts counts_{};
    };

}// namespace transmission_nets::core::samplers

#endif//TRANSMISSION_NETWORKS_APP_ABSTRACTSAMPLER_H
