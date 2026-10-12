//
// The console progress adapter. See Progress.h.
//

#include "Progress.h"

#include <fmt/format.h>

#include <cmath>
#include <optional>
#include <ostream>

namespace transmission_nets::impl::Model {
    namespace {
        std::string phaseName(const Phase phase) {
            return phase == Phase::BurnIn ? "burn-in" : "sampling";
        }

        // "burn-in 150/2000": the Phase, and the Sampling iterations run out of the Run's total.
        std::string where(const Phase phase, const RunPosition& at) {
            return fmt::format("{} {}/{}", phaseName(phase), at.iteration, at.totalIterations());
        }

        // Sampling iterations per second, if any time was spent.
        std::optional<double> rate(const std::int64_t iterations, const double seconds) {
            return seconds > 0 ? std::optional(static_cast<double>(iterations) / seconds) : std::nullopt;
        }

        std::string rateText(const std::optional<double> perSecond) {
            return perSecond ? fmt::format("{:.2f}", *perSecond) : "-";
        }

        // The time left at `perSecond`: "45s", "2m28s" or "2h06m".
        std::string eta(const RunPosition& at, const std::optional<double> perSecond) {
            const auto remaining = at.remaining();
            if (remaining == 0) {
                return "0s";
            }
            if (!perSecond or *perSecond <= 0) {
                return "?";
            }
            const auto s = std::llround(static_cast<double>(remaining) / *perSecond);
            if (s < 60) {
                return fmt::format("{}s", s);
            }
            if (s < 3600) {
                return fmt::format("{}m{:02}s", s / 60, s % 60);
            }
            return fmt::format("{}h{:02}m", s / 3600, s % 3600 / 60);
        }

        std::string logDensities(const double prior, const double likelihood) {
            return fmt::format("prior/likelihood {:.2f} / {:.2f}", prior, likelihood);
        }
    }// namespace

    ConsoleProgress::ConsoleProgress(std::ostream& out) : out_(out) {}

    void ConsoleProgress::line(const std::string& text) {
        out_ << text << '\n';
        out_.flush();
    }

    void ConsoleProgress::notice(const std::string& text) {
        line(text);
    }

    void ConsoleProgress::started(const RunPosition& at, const bool resumed) {
        const auto total = at.totalIterations();
        if (resumed) {
            line(fmt::format("Resuming at iteration {}/{} (Window {} committed).", at.iteration, total, at.lastClosedWindow));
        } else {
            line(fmt::format("Starting at iteration {}/{}.", at.iteration, total));
        }
    }

    void ConsoleProgress::windowClosed(const RunPosition& at, const WindowSummary& window) {
        samplingIterations_ += window.samplingIterations;
        seconds_ += window.seconds;
        const auto windowRate = rate(window.samplingIterations, window.seconds);
        auto text             = fmt::format("{}  Window {}  ETA {}  {} iterations/s (Run {})  {}", where(window.phase, at), window.window, eta(at, windowRate),
                                            rateText(windowRate), rateText(rate(samplingIterations_, seconds_)), logDensities(window.hotPrior, window.hotLikelihood));
        if (!window.swaps.empty()) {
            text += "  swap acceptance";
            for (const auto& [accepted, attempted] : window.swaps) {
                text += attempted == 0 ? std::string(" -") : fmt::format(" {:.2f}", static_cast<double>(accepted) / static_cast<double>(attempted));
            }
        }
        if (adaptations_ > 0) {
            text += fmt::format("  ladder {:.3g} (adapted {}\u00d7)", fmt::join(ladder_, ", "), adaptations_);
            adaptations_ = 0;
        }
        line(text);
    }

    // Adaptations come every few burn-in Sampling iterations, too often for a line each: the next Window line sums them up.
    void ConsoleProgress::ladderAdapted(const RunPosition&, const TemperatureAdaptation& adaptation) {
        if (adaptation.applied) {
            ladder_ = adaptation.temperatures;
            ++adaptations_;
        }
    }

    void ConsoleProgress::heartbeat(const RunPosition& at, const Heartbeat& beat) {
        line(fmt::format("{}  heartbeat  ETA {}  {}", where(at.phase, at), eta(at, rate(beat.samplingIterations, beat.seconds)), logDensities(beat.hotPrior, beat.hotLikelihood)));
    }

    void ConsoleProgress::stopped(const RunPosition& at, const std::string&) {
        line(fmt::format("Stopped at iteration {} (Window {} committed).", at.iteration, at.lastClosedWindow));
    }

    void ConsoleProgress::ended(const RunPosition& at) {
        line(fmt::format("Completed at iteration {} (Window {} committed).", at.iteration, at.lastClosedWindow));
    }
}// namespace transmission_nets::impl::Model
