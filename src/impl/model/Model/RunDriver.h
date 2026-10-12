//
// The Run driver: the one place that advances a Run from where it stands to its end or to a stop. One per Run. It runs
// burn-in and sampling through one iteration body, moves the Run from burn-in to sampling with the final adaptation,
// forwards ladder adaptations, takes the Run's Checkpoints, and decides when a stop request takes effect. It reports
// the Run's progress, and is the one place that decides when a heartbeat is due. It is not stepwise: a caller makes one
// call, and the Run driver returns once the Run has completed or stopped.
//
// Checkpoints. A complete Window commits once the next Sampling iteration has been recorded, so the Run driver
// captures the Checkpoint before running that iteration: the Checkpoint holds the State the Window ends at. A Window
// the Run's end commits ends at the State as it is then. Each Checkpoint is written to <store>/checkpoint.cbor when its
// Window commits, and the one whose Window holds burn-in end to <store>/checkpoint-burn-in.cbor as well.
//

#ifndef TRANSMISSION_NETWORKS_APP_RUNDRIVER_H
#define TRANSMISSION_NETWORKS_APP_RUNDRIVER_H

#include "Checkpoint.h"
#include "ParentSetIds.h"
#include "Progress.h"
#include "RunPosition.h"
#include "RunRecorder.h"

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <mutex>
#include <optional>
#include <string>
#include <variant>

namespace transmission_nets::impl::Model {
    // How long a Run is, in Sampling iterations per Phase.
    struct RunPlan {
        std::int64_t burnIn  = 0;
        std::int64_t samples = 0;
    };

    // Where a Run starts: a fresh Run at its starting State, or a Run resumed from its Checkpoint, where the
    // Checkpoint left it.
    struct FreshRun {};
    struct ResumedRun {
        RunPosition position;
        ParentSetIdDictionary parentSetIds{};
    };
    using RunStart = std::variant<FreshRun, ResumedRun>;

    // How a Run the Run driver advanced ended: it ran to the end of its plan, or a stop request ended it, naming the
    // first signal.
    struct Completed {};
    struct Stopped {
        std::string signal;
    };
    using RunOutcome = std::variant<Completed, Stopped>;

    class RunDriver {
    public:
        // Polled between Sampling iterations: the first signal once the Run has been asked to stop, none before.
        using StopQuery = std::function<std::optional<std::string>()>;

        // Advances `ensemble`, whose top-level engine is `runRng`, records it through `recorder`, writes its
        // Checkpoints, stamped with `identity`, into `store`, and reports its progress to `progress`. `clock` decides
        // when a heartbeat is due; give it the recorder's.
        RunDriver(Ensemble& ensemble, const checkpoint::Engine& runRng, RunRecorder& recorder, std::filesystem::path store, checkpoint::Identity identity,
                  Progress& progress, RunRecorder::Clock clock = std::chrono::system_clock::now);

        // Runs the Run from `start` to the end of `plan`, or until `stop` asks it to stop before a Sampling iteration.
        // Records the Run's start or resumption, burn-in end, every Sampling iteration, any first interrupt and the
        // Run's end. Call once.
        RunOutcome drive(const RunStart& start, const RunPlan& plan, const StopQuery& stop);

        // For a second interrupt, from any thread: ends the Run in the recorder, recording the first interrupt too
        // unless drive() already has, and keeps the recorder from then on, so nothing is recorded after the Run has
        // ended. Returns with the Run driver blocked for good; the caller ends the process. Records nothing outside
        // drive().
        void forceEnd(const std::string& firstSignal, const std::string& signal);

    private:
        Ensemble& ensemble_;
        const checkpoint::Engine& runRng_;
        RunRecorder& recorder_;
        std::filesystem::path store_;
        checkpoint::Identity identity_;
        Progress& progress_;
        RunRecorder::Clock clock_;

        // Guards the recorder between drive() and forceEnd().
        std::mutex lock_;
        bool driving_       = false;
        bool firstRecorded_ = false;
    };
}// namespace transmission_nets::impl::Model

#endif//TRANSMISSION_NETWORKS_APP_RUNDRIVER_H
