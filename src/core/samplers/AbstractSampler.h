//
// Created by Maxwell Murphy on 3/25/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ABSTRACTSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_ABSTRACTSAMPLER_H

#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "core/computation/PartialLikelihood.h"
#include "core/samplers/SamplerEnum.h"

namespace transmission_nets::core::samplers {
    using Likelihood = core::computation::Likelihood;

    // Metropolis-Hastings proposals a sampler has made and accepted since it was built (or since its counts were
    // restored), over every Sampling iteration of both Phases.
    struct SamplerCounts {
        std::uint64_t accepted = 0;
        std::uint64_t proposed = 0;

        bool operator==(const SamplerCounts&) const = default;
    };

    // One proposal scale and the proposals made at it, which the adaptation rule reads.
    struct ProposalScale {
        double variance = 0;
        SamplerCounts counts{};

        bool operator==(const ProposalScale&) const = default;
    };

    // Everything a sampler's adaptation rule reads: its proposal scales (one for a random walk, one per Allele
    // frequency for SALT, none for a sampler that doesn't adapt) and its update count, which sets the step size.
    struct AdaptationState {
        std::vector<ProposalScale> scales{};
        std::uint64_t updates = 0;

        bool operator==(const AdaptationState&) const = default;
    };

    class AbstractSampler {
    public:
        virtual ~AbstractSampler()     = default;
        virtual void update() noexcept = 0;
        virtual void adapt() noexcept {};
        virtual void adapt([[maybe_unused]] unsigned int idx) noexcept {};

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

        // Call once per proposal, after the accept/reject decision.
        void countProposal(const bool accepted) noexcept {
            ++counts_.proposed;
            counts_.accepted += accepted;
        }

    private:
        SamplerCounts counts_{};
    };



}// namespace transmission_nets::core::samplers


#endif//TRANSMISSION_NETWORKS_APP_ABSTRACTSAMPLER_H
