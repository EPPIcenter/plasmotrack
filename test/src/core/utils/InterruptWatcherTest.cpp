//
// The interrupt watcher: a first interrupt asks the Run to stop at its next Sampling iteration; a second one runs the
// force action (the CLI flushes the Run recorder there) and exits. SIGUSR1 and SIGUSR2 never force.
//
// Each case runs in a death-test child, because the watcher blocks the interrupt signals in the process for good.
//

#include "gtest/gtest.h"

#include "core/utils/InterruptWatcher.h"

#include <signal.h>
#include <unistd.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <thread>

using transmission_nets::core::utils::InterruptWatcher;

namespace {
    void waitUntilInterrupted(const InterruptWatcher& watcher) {
        while (!watcher.interrupted()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
    }

    [[noreturn]] void idle() {
        for (;;) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
        }
    }
}// namespace

TEST(InterruptWatcherDeathTest, FirstInterruptOnlyMarksTheRunInterrupted) {
    EXPECT_EXIT(
            {
                InterruptWatcher watcher([](const std::string&) { std::fputs("forced\n", stderr); });
                ::kill(::getpid(), SIGINT);
                waitUntilInterrupted(watcher);
                std::this_thread::sleep_for(std::chrono::milliseconds(50));
                std::fputs("still running\n", stderr);
                std::exit(watcher.firstSignalNumber() == SIGINT ? 7 : 8);
            },
            ::testing::ExitedWithCode(7), "^still running\n$");
}

TEST(InterruptWatcherDeathTest, SecondInterruptRunsTheForceActionThenExits) {
    EXPECT_EXIT(
            {
                InterruptWatcher watcher([](const std::string&) { std::fputs("forced\n", stderr); });
                ::kill(::getpid(), SIGINT);
                waitUntilInterrupted(watcher);
                ::kill(::getpid(), SIGTERM);
                idle();
            },
            ::testing::ExitedWithCode(128 + SIGTERM), "forced");
}

TEST(InterruptWatcherDeathTest, SecondSigintExitsWith130AfterTheForceAction) {
    EXPECT_EXIT(
            {
                InterruptWatcher watcher([](const std::string&) { std::fputs("forced\n", stderr); });
                ::kill(::getpid(), SIGINT);
                waitUntilInterrupted(watcher);
                ::kill(::getpid(), SIGINT);
                idle();
            },
            ::testing::ExitedWithCode(130), "^forced\n$");
}

TEST(InterruptWatcherDeathTest, RepeatedUserSignalsNeverForce) {
    EXPECT_EXIT(
            {
                InterruptWatcher watcher([](const std::string&) { std::fputs("forced\n", stderr); });
                ::kill(::getpid(), SIGUSR1);
                waitUntilInterrupted(watcher);
                ::kill(::getpid(), SIGUSR2);
                std::this_thread::sleep_for(std::chrono::milliseconds(50));
                std::fputs("still running\n", stderr);
                std::exit(7);
            },
            ::testing::ExitedWithCode(7), "^still running\n$");
}

TEST(InterruptWatcherDeathTest, StopsWithoutAnInterrupt) {
    // Constructed and destroyed in a child too: the destructor must join its thread rather than hang.
    EXPECT_EXIT(
            {
                { InterruptWatcher watcher([](const std::string&) {}); }
                std::exit(5);
            },
            ::testing::ExitedWithCode(5), "");
}
