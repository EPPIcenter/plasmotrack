//
// Turns interrupt signals into Run decisions without running anything inside a signal handler. The constructor blocks
// SIGINT, SIGQUIT, SIGABRT, SIGTERM, SIGUSR1 and SIGUSR2 in the calling thread, so every thread started after it
// inherits the block, and a thread of its own takes them with sigwait. Its reactions therefore run in ordinary thread
// context, where I/O, fmt and gzip are safe.
//
// The first signal marks the Run interrupted. A later SIGINT, SIGQUIT, SIGABRT or SIGTERM runs the force action and
// ends the process with exit code 128 + the signal number (130 for SIGINT); SIGUSR1 and SIGUSR2 never force. The force action runs while other threads keep
// running, so anything it touches must be guarded against them.
//
// Construct it at the start of main, before any other thread exists.
//

#ifndef TRANSMISSION_NETWORKS_APP_INTERRUPTWATCHER_H
#define TRANSMISSION_NETWORKS_APP_INTERRUPTWATCHER_H

#include <atomic>
#include <functional>
#include <string>
#include <thread>

namespace transmission_nets::core::utils {
    class InterruptWatcher {
    public:
        // `onForce` is handed the name of the signal that forces, such as "SIGTERM".
        explicit InterruptWatcher(std::function<void(const std::string& signal)> onForce);
        ~InterruptWatcher();
        InterruptWatcher(const InterruptWatcher&)            = delete;
        InterruptWatcher& operator=(const InterruptWatcher&) = delete;

        // Whether an interrupt has been received.
        [[nodiscard]] bool interrupted() const;

        // The name of the first interrupt signal, such as "SIGINT". Empty until one is received.
        [[nodiscard]] std::string firstSignal() const;

        // The number of the first interrupt signal, such as SIGINT; 0 until one is received. A Run the first interrupt
        // stopped exits with 128 + this number once it has shut down cleanly.
        [[nodiscard]] int firstSignalNumber() const;

    private:
        void watch();

        std::function<void(const std::string&)> onForce_;
        std::atomic<const char*> firstSignal_{nullptr};
        std::atomic<int> firstSignalNumber_{0};
        std::atomic<bool> interrupted_{false};
        std::atomic<bool> stopping_{false};
        std::thread thread_;
    };
}// namespace transmission_nets::core::utils

#endif//TRANSMISSION_NETWORKS_APP_INTERRUPTWATCHER_H
