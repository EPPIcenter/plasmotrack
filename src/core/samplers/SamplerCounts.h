//
// What a sampler reports about its proposals: how many it made and accepted, and the state its adaptation rule reads.
//

#ifndef TRANSMISSION_NETWORKS_APP_SAMPLERCOUNTS_H
#define TRANSMISSION_NETWORKS_APP_SAMPLERCOUNTS_H

#include <cstdint>
#include <vector>

namespace transmission_nets::core::samplers {

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

}// namespace transmission_nets::core::samplers

#endif//TRANSMISSION_NETWORKS_APP_SAMPLERCOUNTS_H
