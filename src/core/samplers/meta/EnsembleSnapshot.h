//
// What replica exchange hands over about its ensemble: what occupies each temperature slot, how often neighbouring
// slots swapped, and what its temperature adaptation decided. The Run recorder takes these as they are.
//

#ifndef TRANSMISSION_NETWORKS_APP_ENSEMBLESNAPSHOT_H
#define TRANSMISSION_NETWORKS_APP_ENSEMBLESNAPSHOT_H

#include "core/samplers/AbstractSampler.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace transmission_nets::core::samplers {

    // What occupies one temperature slot at the time of a call. Slot 0 is the Hot chain; temperatures fall with the
    // slot.
    struct TemperatureSlot {
        std::size_t replica;
        double temperature;
        double prior;
        double likelihood;
    };

    // Swaps between temperature slots `pair` and `pair + 1`, counted over the whole Run in both Phases.
    struct SwapCounts {
        std::uint64_t accepted  = 0;
        std::uint64_t attempted = 0;

        bool operator==(const SwapCounts&) const = default;
    };

    // What a temperature adaptation decided. `temperatures` is the ladder by temperature slot after it: the new one
    // when `applied`, the unchanged one otherwise. `barriers` holds, per neighbouring pair of slots, the mean swap
    // rejection probability since the previous adaptation, which set the new ladder.
    struct TemperatureAdaptation {
        bool applied;
        std::vector<double> temperatures;
        std::vector<double> barriers;
    };

    // Where adaptation ended when burn-in ends.
    struct FinalAdaptation {
        std::vector<AdaptationState> samplers;// the Hot chain's, one per sampler id in registration order
        std::vector<double> temperatures;     // the ladder by temperature slot
    };

    // The replica running at the target posterior, at the time of a call.
    template<typename State, typename Model>
    struct HotChain {
        State& state;
        Model& model;
    };

    // The ensemble after a Sampling iteration: the Hot chain, every temperature slot, and cumulative counts. Counts
    // never fall between snapshots of one Run.
    template<typename State, typename Model>
    struct EnsembleSnapshot {
        HotChain<State, Model> hot;
        std::vector<TemperatureSlot> slots;                // one per chain, by temperature slot
        std::vector<std::vector<SamplerCounts>> acceptance;// by replica, then by sampler id in registration order
        std::vector<SwapCounts> swaps;                     // one per neighbouring pair of temperature slots
    };
}// namespace transmission_nets::core::samplers

#endif//TRANSMISSION_NETWORKS_APP_ENSEMBLESNAPSHOT_H
