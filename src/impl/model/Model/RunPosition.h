//
// Where a Run stands between Sampling iterations: what the Run recorder reports when it commits a Window, what a
// checkpoint records, and where a resumed Run picks up.
//

#ifndef TRANSMISSION_NETWORKS_APP_RUNPOSITION_H
#define TRANSMISSION_NETWORKS_APP_RUNPOSITION_H

#include <cstdint>

namespace transmission_nets::impl::Model {

    // Which part of a Run a Sampling iteration belongs to.
    enum class Phase { BurnIn, Sampling };

    struct RunPosition {
        std::int64_t iteration         = 0;// the last Sampling iteration recorded, counted within the Run; 0 is the starting State
        Phase phase                    = Phase::BurnIn;// the Phase of the next Sampling iteration: sampling once burn-in end is recorded
        std::int64_t lastClosedWindow  = 0;// the newest committed Window
        std::int64_t burnInRemaining   = 0;
        std::int64_t samplingRemaining = 0;

        // Sampling iterations still to run, in both Phases.
        [[nodiscard]] std::int64_t remaining() const { return burnInRemaining + samplingRemaining; }
        // Sampling iterations the Run plans in all, run or not.
        [[nodiscard]] std::int64_t totalIterations() const { return iteration + remaining(); }

        bool operator==(const RunPosition&) const = default;
    };
}// namespace transmission_nets::impl::Model

#endif//TRANSMISSION_NETWORKS_APP_RUNPOSITION_H
