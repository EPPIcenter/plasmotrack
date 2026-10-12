//
// How the Run driver reports a Run's progress: when the Run starts, each Window it closes, each ladder adaptation, a
// heartbeat when nothing has been reported for a while, and how the Run ends. The console adapter writes these as
// lines; tests record them. Every call says where the Run stands.
//

#ifndef TRANSMISSION_NETWORKS_APP_PROGRESS_H
#define TRANSMISSION_NETWORKS_APP_PROGRESS_H

#include "RunPosition.h"
#include "RunRecorder.h"

#include <chrono>
#include <cstdint>
#include <iosfwd>
#include <optional>
#include <string>
#include <vector>

namespace transmission_nets::impl::Model {

    // What the Run driver knows between Sampling iterations when it sends a heartbeat. The counts are this process's:
    // since the Run started or resumed.
    struct Heartbeat {
        double hotPrior      = 0;
        double hotLikelihood = 0;
        std::int64_t samplingIterations = 0;
        double seconds                  = 0;
    };

    class Progress {
    public:
        virtual ~Progress() = default;

        // The Run starts at `at`, or resumes there.
        virtual void started(const RunPosition& at, bool resumed) = 0;
        // A Window committed; `at` is where the Run stands once it has.
        virtual void windowClosed(const RunPosition& at, const WindowSummary& window) = 0;
        // Replica exchange adapted its temperature ladder after the Sampling iteration `at` ends.
        virtual void ladderAdapted(const RunPosition& at, const TemperatureAdaptation& adaptation) = 0;
        // Nothing but ladder adaptations has been reported for the heartbeat interval.
        virtual void heartbeat(const RunPosition& at, const Heartbeat& beat) = 0;
        // A stop request ended the Run at `at`, for `signal`.
        virtual void stopped(const RunPosition& at, const std::string& signal) = 0;
        // The Run ran to the end of its plan.
        virtual void ended(const RunPosition& at) = 0;
    };

    // How long the Run driver lets pass without a report before it sends a heartbeat.
    inline constexpr std::chrono::seconds HEARTBEAT_INTERVAL{60};

    // The console: whole lines, never redrawn, flushed at once so a Run sent to a file shows them. One line per call,
    // except ladder adaptations, which the next Window line sums up.
    class ConsoleProgress final : public Progress {
    public:
        explicit ConsoleProgress(std::ostream& out);

        void started(const RunPosition& at, bool resumed) override;
        void windowClosed(const RunPosition& at, const WindowSummary& window) override;
        void ladderAdapted(const RunPosition& at, const TemperatureAdaptation& adaptation) override;
        void heartbeat(const RunPosition& at, const Heartbeat& beat) override;
        void stopped(const RunPosition& at, const std::string& signal) override;
        void ended(const RunPosition& at) override;

        // A line before the Run starts: what the Run recipe notes, the null model, the seed drawn.
        void notice(const std::string& text);

    private:
        void line(const std::string& text);

        std::ostream& out_;
        // Over the Windows this process has closed: the Run's rate.
        std::int64_t samplingIterations_ = 0;
        double seconds_                  = 0;
        // The ladder as last adapted, and how many adaptations changed it since the last Window line.
        std::vector<double> ladder_{};
        std::int64_t adaptations_ = 0;
    };
}// namespace transmission_nets::impl::Model

#endif//TRANSMISSION_NETWORKS_APP_PROGRESS_H
