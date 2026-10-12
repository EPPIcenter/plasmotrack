//
// The Run recorder writing a store: the one place that decides what of a Run reaches its output and when. One per Run,
// never one per chain. It writes Trace rows on the `-t` grid in both Phases, Blocks of acceptance and swap counts over
// every Sampling iteration, and the Run events the caller reports. It groups them into Windows counted in Sampling
// iterations and hands each Window to a StoreSink to commit. Rows go to the sink as they are written, and a row group
// may end after each Trace row, so the sink can bound what a Window holds in memory. It reads nothing from replica exchange
// or the samplers itself: the caller hands over what they expose.
//
// Every committed Window covers at least one Sampling iteration. A Window that fills, or that burn-in end is recorded
// in, stays open until the next Sampling iteration is recorded or the Run ends, so Run events raised in between go
// into it.
//

#ifndef TRANSMISSION_NETWORKS_APP_RUNRECORDER_H
#define TRANSMISSION_NETWORKS_APP_RUNRECORDER_H

#include "Model.h"
#include "ParentSetIds.h"
#include "Refusal.h"
#include "RunPosition.h"
#include "State.h"

#include "core/io/StoreSink.h"
#include "core/samplers/AbstractSampler.h"
#include "core/samplers/meta/EnsembleSnapshot.h"

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace transmission_nets::impl::Model {

    // How the recorder thins and groups rows. Recorded in the Run manifest.
    struct RecorderOptions {
        std::int64_t thin     = 10;
        std::int64_t window   = 100;// in Sampling iterations, a multiple of `thin`
        double parentSetFloor = 1e-9;

        bool operator==(const RecorderOptions&) const = default;
    };

    // The names of a Run's inputs, as the Run manifest keys their hashes.
    namespace input_names {
        inline constexpr const char* INPUT            = "input";
        inline constexpr const char* SYMPTOMATIC_IDP  = "symptomatic_idp";
        inline constexpr const char* ASYMPTOMATIC_IDP = "asymptomatic_idp";
    }// namespace input_names

    // The identity and recipe of a Run, as the CLI knows it. The recorder adds its own options when it writes the
    // Run manifest.
    struct RunManifest {
        std::string runId;
        std::string startedAt;
        long seed = 0;
        std::vector<std::string> cliArguments{};
        std::string engineVersion;
        std::string gitCommit;
        std::map<std::string, std::string> inputSha256{};// by input_names
        int chains = 1;
        std::optional<double> gradient{};// `-g`; none for a Run started from a checkpoint, whose ladder came from there
        std::int64_t burnIn  = 0;
        std::int64_t samples = 0;
        bool nullModel       = false;
        std::optional<nlohmann::json> startedFrom{};

        bool operator==(const RunManifest&) const = default;
    };

    // What the store of a Run records about it: its Run manifest and the recorder's options.
    struct StoredRun {
        RunManifest manifest;
        RecorderOptions options;

        bool operator==(const StoredRun&) const = default;
    };

    // The Run recorder owns the Run manifest's JSON both ways: this reads what it writes. Throws nlohmann::json
    // exceptions for JSON that lacks a field.
    StoredRun storedRunFromJson(const nlohmann::json& j);

    // The Run the store at `store` holds, or none if it holds no Run manifest. Throws Refused for a Run manifest it
    // can't read.
    std::optional<StoredRun> readStoredRun(const std::filesystem::path& store);

    // The ensemble's view, as replica exchange hands it over; see core/samplers/meta/EnsembleSnapshot.h.
    using core::samplers::FinalAdaptation;
    using core::samplers::SwapCounts;
    using core::samplers::TemperatureAdaptation;
    using core::samplers::TemperatureSlot;
    using HotChain         = core::samplers::HotChain<State, Model>;
    using EnsembleSnapshot = core::samplers::EnsembleSnapshot<State, Model>;

    // A Run resumed from its checkpoint: where the checkpoint left it, and what the resumed Run plans now.
    struct Resumption {
        RunPosition position;
        std::string checkpoint;// the checkpoint's path
        std::int64_t samples;  // sampling iterations planned now, which may have grown
        ParentSetIdDictionary parentSetIds{};
    };

    // What a committed Window held, for reporting progress. Its time and Sampling iterations are those of its
    // window_closed Run event, and its swap counts its swaps Block.
    struct WindowSummary {
        std::int64_t window = 0;
        Phase phase         = Phase::BurnIn;
        std::int64_t samplingIterations = 0;
        double seconds                  = 0;
        // The Hot chain's log prior and log likelihood at the Window's last Sampling iteration.
        double hotPrior      = 0;
        double hotLikelihood = 0;
        // Per neighbouring pair of temperature slots, the swaps over the Window.
        std::vector<SwapCounts> swaps{};

        bool operator==(const WindowSummary&) const = default;
    };

    // A committed Window: where the Run stands once its marker is written, which is what a checkpoint taken then
    // covers, whether burn-in end was recorded in it, and what it held.
    struct WindowCommit {
        RunPosition position;
        const ParentSetIdDictionary& parentSetIds;
        bool burnInEnded;
        WindowSummary summary{};
    };

    // Which interrupt of the Run a signal is.
    enum class Interrupt { First, Second };

    class RunRecorder {
    public:
        // Writes the Run manifest and the structure tables. `samplerIds` are the `<type>/<target>` ids of the samplers
        // whose counts each snapshot carries, in the same order. `options.window` must be a positive multiple of
        // `options.thin`, as the Run recipe checks.
        //
        // `clock` gives the wall time of each Run event and the time spent in each Window.
        using Clock = std::function<std::chrono::system_clock::time_point()>;
        RunRecorder(const State& structure, std::vector<std::string> samplerIds, std::unique_ptr<core::io::StoreSink> sink, RecorderOptions options, const RunManifest& manifest,
                      Clock clock = std::chrono::system_clock::now);
        RunRecorder(RunRecorder&&) noexcept;
        ~RunRecorder();

        // Called once the Window's marker is written, whenever the Run's position is one a resumed Run can continue
        // from: not after the Window a second interrupt closes, since the State may then be mid Sampling iteration.
        // A Window that fills commits when the next Sampling iteration is recorded, after that iteration has run, so
        // a checkpoint taken here must hold the State as it was when windowComplete() last became true; see
        // windowComplete().
        using WindowCommitted = std::function<void(const WindowCommit&)>;
        void onWindowCommitted(WindowCommitted committed);

        // The starting state: Sampling iteration 0, burn-in, one row on every Trace. Call once, before any iteration.
        void recordStart(const EnsembleSnapshot& ensemble);

        // In place of recordStart, for a Run resumed from its checkpoint: first discards every Window newer than the
        // checkpoint's, which a Run killed between a Window's marker and its checkpoint leaves behind. Iteration and
        // Window numbering continue after `resumption.position`, and Blocks count from the counts in `ensemble`. Writes
        // no Trace row. Call once, before any iteration.
        void recordResume(const Resumption& resumption, const EnsembleSnapshot& ensemble);

        // A finished Sampling iteration in the current Phase. First commits the in-flight Window if it is complete,
        // then writes Trace rows on the `-t` grid. The Window is complete once it holds `window` Sampling iterations.
        // Returns whether Trace rows were written.
        bool recordIteration(const EnsembleSnapshot& ensemble);

        // Replica exchange adapted its temperature ladder after the last Sampling iteration recorded.
        void recordTemperatureAdaptation(const TemperatureAdaptation& adaptation);

        // Burn-in is over and adaptation stops at `adaptation`: the in-flight Window is complete, and the sampling
        // Phase starts with the next one. Throws std::invalid_argument unless `adaptation` has one state per sampler
        // id, and std::logic_error if no Sampling iteration has been recorded since the Run resumed.
        void recordBurnInEnd(const FinalAdaptation& adaptation);

        // A signal asked the Run to stop. After the first, the caller finishes the Sampling iteration in flight and
        // ends the Run with recordEnd. The second ends the Run here: it closes the in-flight Window, and the recorder
        // records nothing more. Safe to call with a Sampling iteration in flight only if no other call is.
        void recordInterrupt(const std::string& signal, Interrupt which);

        // The Run is over, completed or after an interrupt: commits the in-flight Window. A resumed Run that recorded
        // no Sampling iteration commits nothing, so the store stays as its checkpoint left it. Does nothing once a
        // second interrupt has ended the Run.
        void recordEnd();

        // Whether the in-flight Window is complete: full, or holding burn-in end. It commits when the next Sampling
        // iteration is recorded, so a caller that takes checkpoints takes the State now, before running that
        // iteration.
        [[nodiscard]] bool windowComplete() const { return complete_; }

        // The Phase of the last Sampling iteration recorded.
        [[nodiscard]] Phase phase() const { return phase_; }

    private:
        struct WindowRows;

        void writeStructure(const State& structure);
        void writeTraces(const EnsembleSnapshot& ensemble);
        void takeCounts(const EnsembleSnapshot& ensemble);
        std::vector<SwapCounts> writeBlocks();
        void checkFinite(const EnsembleSnapshot& ensemble);
        void raise(const std::string& type, nlohmann::json payload);
        void handOverRows();
        void commitWindow(std::chrono::system_clock::time_point closedAt);
        [[nodiscard]] RunPosition position() const;

        std::vector<std::string> infectionIds_;
        std::vector<std::string> latentParentIds_;
        std::vector<std::string> locusLabels_;
        std::vector<std::string> samplerIds_;
        std::size_t chains_;

        std::unique_ptr<core::io::StoreSink> sink_;
        RecorderOptions options_;
        Clock clock_;
        std::int64_t burnIn_;
        std::int64_t samples_;
        std::chrono::system_clock::time_point startedAt_{};    // when this process began the Run or resumed it
        std::chrono::system_clock::time_point windowOpenedAt_{};
        std::chrono::system_clock::time_point lastRecordedAt_{};// when the last Sampling iteration was recorded
        bool interrupted_ = false;
        bool ended_       = false;
        bool forced_      = false;// a second interrupt ended the Run

        Phase phase_                 = Phase::BurnIn;// of the in-flight Window
        bool burnInEnded_            = false;        // burn-in end is recorded in the in-flight Window
        bool complete_               = false;        // the in-flight Window commits before the next Sampling iteration
        std::int64_t iteration_      = 0;// the last Sampling iteration recorded; 0 is the starting state
        std::int64_t phaseIterations_ = 0;// Sampling iterations recorded in this Phase
        std::int64_t window_         = 0;// the in-flight Window
        std::int64_t windowIterations_ = 0;// Sampling iterations in it, the starting state included
        std::int64_t windowSampled_    = 0;// of them, Sampling iterations after the starting state
        std::unique_ptr<WindowRows> rows_;

        // Cumulative counts as of the last call, and as of the last Block. Acceptance counts only the Hot chain: what
        // each Sampling iteration added to the replica in temperature slot 0 while it ran.
        struct Counts {
            std::vector<core::samplers::SamplerCounts> acceptance;
            std::vector<SwapCounts> swaps;
        };
        Counts latest_{};
        Counts blockStart_{};
        std::vector<std::vector<core::samplers::SamplerCounts>> replicaCounts_;// each replica's counts as of the last call
        std::size_t hotReplica_ = 0;                                            // the replica in slot 0 after the last call
        struct {
            double prior      = 0;
            double likelihood = 0;
        } hot_{};// temperature slot 0's log densities at the last call
        // Per temperature slot, whether its prior and its likelihood were non-finite at the last call.
        std::vector<std::pair<bool, bool>> nonFinite_;
        ParentSetIdDictionary parentSetIds_;
        WindowCommitted windowCommitted_{};
    };
}// namespace transmission_nets::impl::Model

#endif//TRANSMISSION_NETWORKS_APP_RUNRECORDER_H
