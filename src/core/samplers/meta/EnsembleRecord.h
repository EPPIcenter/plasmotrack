//
// What replica exchange carries between Sampling iterations, as a record it makes and restores from itself: the
// temperature slots, ladder, pending swap barriers and swap counts; the top-level engine; and, per replica, its chain's
// State as the caller records it, its engine and its Sampler set. ReplicaExchange::record() makes one, and ReplicaExchange
// has a constructor that restores from one. The Checkpoint is this record's codec.
//

#ifndef TRANSMISSION_NETWORKS_APP_ENSEMBLERECORD_H
#define TRANSMISSION_NETWORKS_APP_ENSEMBLERECORD_H

#include "core/samplers/meta/EnsembleSnapshot.h"
#include "core/samplers/scheduler/Scheduler.h"

#include <cstddef>
#include <string>
#include <vector>

namespace transmission_nets::core::samplers {
    // Everything replica exchange carries between Sampling iterations, apart from each chain's own State, engine and
    // Sampler set.
    struct ExchangeState {
        std::vector<std::size_t> replicas;// replica occupying each temperature slot
        std::vector<double> temperatures; // ladder by temperature slot
        std::vector<double> swapBarriers; // summed swap rejection per pair since the last adaptation
        int barrierIterations  = 0;       // burn-in iterations those sums cover
        int burnInIterations   = 0;       // burn-in iterations run, which decide when the ladder next adapts
        bool secondPairSetNext = false;   // the next iteration tries pairs (1,2),(3,4),... instead of (0,1),...
        std::vector<SwapCounts> swapCounts;

        bool operator==(const ExchangeState&) const = default;
    };

    // The ensemble's record. `ChainState` is a chain's State as the caller records it.
    template<typename ChainState>
    struct EnsembleRecord {
        struct Replica {
            ChainState state;
            std::string rng;          // the chain's engine
            SamplerSetRecord samplers;// its Sampler set: iterations, counts and adaptation state

            bool operator==(const Replica&) const = default;
        };

        std::vector<std::string> samplerIds;// every replica's sampler ids, in registration order
        std::string runRng;                 // the top-level engine, which draws swaps
        std::vector<Replica> replicas;      // by replica, not by temperature slot
        ExchangeState exchange;

        bool operator==(const EnsembleRecord&) const = default;
    };
}// namespace transmission_nets::core::samplers

#endif//TRANSMISSION_NETWORKS_APP_ENSEMBLERECORD_H
