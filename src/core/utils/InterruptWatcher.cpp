//
// The interrupt watcher. See InterruptWatcher.h.
//

#include "InterruptWatcher.h"

#include <fmt/core.h>

#include <pthread.h>
#include <signal.h>

#include <algorithm>
#include <array>
#include <csignal>
#include <cstdio>
#include <cstdlib>

namespace transmission_nets::core::utils {
    namespace {
        struct WatchedSignal {
            int number;
            const char* name;
        };

        // Every signal the watcher takes.
        constexpr std::array<WatchedSignal, 6> WATCHED_SIGNALS{{
                {SIGINT, "SIGINT"},
                {SIGQUIT, "SIGQUIT"},
                {SIGABRT, "SIGABRT"},
                {SIGTERM, "SIGTERM"},
                {SIGUSR1, "SIGUSR1"},
                {SIGUSR2, "SIGUSR2"},
        }};

        sigset_t watchedSignals() {
            sigset_t set;
            sigemptyset(&set);
            for (const auto& [number, _] : WATCHED_SIGNALS) {
                sigaddset(&set, number);
            }
            return set;
        }

        const char* signalName(const int sig) {
            const auto it = std::ranges::find(WATCHED_SIGNALS, sig, &WatchedSignal::number);
            return it != WATCHED_SIGNALS.end() ? it->name : "(unidentified)";
        }
    }// namespace

    InterruptWatcher::InterruptWatcher(std::function<void(const std::string&)> onForce) : onForce_(std::move(onForce)) {
        const auto set = watchedSignals();
        pthread_sigmask(SIG_BLOCK, &set, nullptr);
        // An inherited SIG_IGN (a shell starts background jobs ignoring SIGINT) would discard a signal before sigwait
        // sees it. The default action never runs while the signals are blocked.
        for (const auto& [number, _] : WATCHED_SIGNALS) {
            std::signal(number, SIG_DFL);
        }
        thread_ = std::thread([this] { watch(); });
    }

    InterruptWatcher::~InterruptWatcher() {
        stopping_ = true;
        // Wakes the watcher's sigwait; it sees `stopping_` and returns.
        pthread_kill(thread_.native_handle(), SIGUSR1);
        thread_.join();
    }

    bool InterruptWatcher::interrupted() const {
        return interrupted_;
    }

    std::string InterruptWatcher::firstSignal() const {
        const char* name = firstSignal_;
        return name == nullptr ? "" : name;
    }

    int InterruptWatcher::firstSignalNumber() const {
        return firstSignalNumber_;
    }

    void InterruptWatcher::watch() {
        const auto set = watchedSignals();
        for (;;) {
            // macOS can report 0 for one of two different signals pending at once; it still counts as an interrupt.
            int sig = 0;
            if (sigwait(&set, &sig) != 0 or stopping_) {
                return;
            }
            if (interrupted_ and sig != SIGUSR1 and sig != SIGUSR2) {
                fmt::print("Interrupt signal {} received after an earlier one. Flushing output and exiting...\n", signalName(sig));
                std::fflush(stdout);
                onForce_(signalName(sig));
                std::fflush(stdout);
                std::fflush(stderr);
                // The main thread may be mid-iteration: skip static destructors and atexit handlers it could race.
                // 128 + the signal number is how a shell reports a process the signal ended.
                std::_Exit(128 + sig);
            }
            firstSignal_       = signalName(sig);
            firstSignalNumber_ = sig;
            interrupted_       = true;
            // One whole line, so it never runs into a progress line.
            fmt::print("Interrupt signal {} received. Stopping after the Sampling iteration in flight...\n", signalName(sig));
            std::fflush(stdout);
        }
    }
}// namespace transmission_nets::core::utils
