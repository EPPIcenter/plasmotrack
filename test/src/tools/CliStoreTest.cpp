//
// The CLI's store end to end under real signals: a first interrupt ends the Run cleanly and a second ends it at once.
// Each test runs the CLI binary on the example input and reads the store back with the Arrow C++ reader. How the Run
// driver ends a Run for a stop request or a forced end is tested in-process in RunDriverTest.
//

#include "gtest/gtest.h"

#include "CliRunning.h"

#include "impl/model/Model/Checkpoint.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

using namespace transmission_nets::test;
namespace checkpoint = transmission_nets::impl::Model::checkpoint;

namespace {
    // A long Run: it ends only by a signal within any test's time.
    const std::vector<std::string> LONG_RUN{"--seed", "7", "-b", "0", "-s", "100000000", "-t", "1", "--window", "10"};
}// namespace

// A first interrupt lets the Sampling iteration in flight finish, then ends the Run cleanly: it records the interrupt,
// commits the in-flight Window and its checkpoint, and exits with 128 + the signal, as a shell expects of a process a
// signal stopped.
TEST(CliStore, FirstInterruptEndsTheRunAndClosesTheInFlightWindow) {
    const auto work  = makeWork("interrupt");
    const auto store = work / "out";
    const pid_t pid  = startCli(work, LONG_RUN);
    ASSERT_TRUE(waitForClosedWindows(pid, store, 3)) << "console output in " << work / "cli.txt";
    ::kill(pid, SIGINT);
    EXPECT_EQ(waitFor(pid), 128 + SIGINT) << "console output in " << work / "cli.txt";
    const std::int64_t newest = readJson(store / "_closed" / closedWindows(store).back()).at("window");
    EXPECT_EQ(checkpoint::load(store / "checkpoint.cbor").position.lastClosedWindow, newest);

    expectEveryClosedWindowReadsBack(store);
    const auto e = events(store);
    const auto t = types(e);
    ASSERT_GE(e.size(), 3u);
    EXPECT_EQ(std::vector(t.end() - 3, t.end()), (std::vector<std::string>{"interrupted", "run_ended", "window_closed"}));
    EXPECT_EQ(e[e.size() - 3].payload, (nlohmann::json{{"signal", "SIGINT"}, {"second", false}}));
    EXPECT_EQ(e[e.size() - 2].payload.at("reason"), "interrupted");

    if (!HasFailure()) fs::remove_all(work);
}

// A second interrupt ends the Run at once, mid Sampling iteration: it closes the in-flight Window and exits with
// 128 + the signal. Sixty-four chains make a Sampling iteration long (about 10 ms). SIGINT goes once Window 1's
// checkpoint is written, which is the last thing the Run does before it starts the next Sampling iteration, rather than
// as soon as the Window's marker appears, while the Run is still between iterations. SIGTERM goes once the CLI reports
// the first interrupt, so the CLI takes the two signals one at a time and in order. If the Run still ended after the
// first interrupt (the CLI finished its Sampling iteration before SIGTERM came), the test tries again.
TEST(CliStore, SecondInterruptClosesTheInFlightWindowAndExitsWithTheSignal) {
    auto options = LONG_RUN;
    options.insert(options.end(), {"-n", "64", "-g", "0.1"});
    for (int attempt = 0; attempt < 20; ++attempt) {
        const auto work  = makeWork("second-interrupt");
        const auto store = work / "out";
        const pid_t pid  = startCli(work, options);
        // Window 0 holds the starting State and burn-in end and commits with the first sampling iteration; Window 1
        // with the eleventh.
        ASSERT_TRUE(waitForClosedWindows(pid, store, 2)) << "console output in " << work / "cli.txt";
        const std::int64_t committed = readJson(store / "_closed" / closedWindows(store).back()).at("window");
        ASSERT_TRUE(waitUntil(pid, [&] { return checkpoint::load(store / "checkpoint.cbor").position.lastClosedWindow >= committed; }))
                << "console output in " << work / "cli.txt";
        ::kill(pid, SIGINT);
        // False only if the CLI has exited, which the events below show.
        waitForConsoleOutput(pid, work, "Interrupt signal SIGINT received");
        ::kill(pid, SIGTERM);
        const int exitCode = waitFor(pid);
        const auto e       = events(store);
        const auto t       = types(e);
        if (std::ranges::count(t, "interrupted") < 2) {
            continue;// the Run had ended after the first interrupt when SIGTERM came
        }
        EXPECT_EQ(exitCode, 128 + SIGTERM) << "console output in " << work / "cli.txt";

        expectEveryClosedWindowReadsBack(store);
        ASSERT_GE(e.size(), 4u);
        EXPECT_EQ(std::vector(t.end() - 4, t.end()), (std::vector<std::string>{"interrupted", "interrupted", "run_ended", "window_closed"}));
        EXPECT_EQ(e[e.size() - 4].payload, (nlohmann::json{{"signal", "SIGINT"}, {"second", false}}));
        EXPECT_EQ(e[e.size() - 3].payload, (nlohmann::json{{"signal", "SIGTERM"}, {"second", true}}));
        // The Window the second interrupt closed has no checkpoint: the State may have been mid Sampling iteration.
        const std::int64_t newest = readJson(store / "_closed" / closedWindows(store).back()).at("window");
        EXPECT_LT(checkpoint::load(store / "checkpoint.cbor").position.lastClosedWindow, newest);
        if (!HasFailure()) fs::remove_all(work);
        return;
    }
    FAIL() << "every Run ended after the first interrupt, before the second came";
}
