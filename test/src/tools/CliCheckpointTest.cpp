//
// Checkpoints through the CLI end to end: --resume continues an interrupted Run from its checkpoint, --start-from
// begins a new Run from one's State, and a checkpoint that can't be written ends the Run. Each test runs the CLI binary
// on the example input and reads the store back with the Arrow C++ reader. When Checkpoints are taken and where a
// resumed Run picks up are tested in-process in RunDriverTest; which checkpoints a Run may resume or start from, in
// RunRecipeTest.
//

#include "gtest/gtest.h"

#include "CliRunning.h"

#include "impl/model/Model/Checkpoint.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using namespace transmission_nets::test;
namespace checkpoint = transmission_nets::impl::Model::checkpoint;

namespace {
    using Phase = checkpoint::Phase;

    // Each closed Window's number, from its marker, in order.
    std::vector<std::int64_t> windowNumbers(const fs::path& store) {
        std::vector<std::int64_t> out;
        for (const auto& w : closedWindows(store)) {
            out.push_back(readJson(store / "_closed" / w).at("window"));
        }
        return out;
    }

    std::vector<std::int64_t> zeroTo(const std::int64_t last, const std::int64_t step = 1) {
        std::vector<std::int64_t> out;
        for (std::int64_t k = 0; k <= last; k += step) {
            out.push_back(k);
        }
        return out;
    }

    std::vector<std::int64_t> traceIterations(const fs::path& store) {
        return readStoreTable(store, "run_scalars").column<std::int64_t>("iteration");
    }

    // A run_scalars column at one iteration.
    double scalarAt(const fs::path& store, const std::string& column, const std::int64_t iteration) {
        const auto t          = readStoreTable(store, "run_scalars");
        const auto iterations = t.column<std::int64_t>("iteration");
        const auto at         = std::ranges::find(iterations, iteration);
        EXPECT_NE(at, iterations.end()) << "no row at iteration " << iteration;
        return t.column<double>(column)[at - iterations.begin()];
    }

    // Resumes the Run in <work>/out with `options`, writing console output to <work>/cli.txt.
    pid_t startResume(const fs::path& work, std::vector<std::string> options) {
        options.insert(options.begin(), {"--resume", (work / "out").string()});
        return spawnCli(work, cliArguments({}, options));
    }

    int resume(const fs::path& work, const std::vector<std::string>& options) {
        return waitFor(startResume(work, options));
    }

    std::vector<Event> eventsOfType(const fs::path& store, const std::string& type) {
        std::vector<Event> out;
        for (const auto& e : events(store)) {
            if (e.type == type) {
                out.push_back(e);
            }
        }
        return out;
    }

    // A long Run: it ends only by a signal within any test's time.
    const std::vector<std::string> LONG_RUN{"--seed", "7", "-b", "0", "-s", "100000000", "-t", "1", "--window", "10"};
    // A short Run: Window 0 (iterations 0-10, with burn-in end), then 1 and 2 (sampling, the Run's end in 2).
    const std::vector<std::string> SHORT_RUN{"--seed", "3", "-b", "10", "-s", "20", "-t", "5", "--window", "10"};
}// namespace

// Interrupted, resumed and interrupted again: one Run, whose iterations and Windows run on without a gap or a repeat.
TEST(CliCheckpoint, AnInterruptedRunResumesWhereItsCheckpointLeftOff) {
    const auto work  = makeWork("resume");
    const auto store = work / "out";
    pid_t pid        = startCli(work, LONG_RUN);
    ASSERT_TRUE(waitForClosedWindows(pid, store, 3)) << consoleOutput(work);
    ::kill(pid, SIGINT);
    ASSERT_EQ(waitFor(pid), 128 + SIGINT) << consoleOutput(work);
    const auto manifest     = readJson(store / "manifest.json");
    const auto interruptedAt = traceIterations(store).back();
    const auto windowsBefore = closedWindows(store).size();

    pid = startResume(work, {});
    ASSERT_TRUE(waitForClosedWindows(pid, store, windowsBefore + 3)) << consoleOutput(work);
    ::kill(pid, SIGINT);
    ASSERT_EQ(waitFor(pid), 128 + SIGINT) << consoleOutput(work);

    expectEveryClosedWindowReadsBack(store);
    const auto windows = windowNumbers(store);
    EXPECT_EQ(windows, zeroTo(windows.back()));
    const auto iterations = traceIterations(store);
    EXPECT_EQ(iterations, zeroTo(iterations.back()));
    EXPECT_GT(iterations.back(), interruptedAt);

    const auto resumed = eventsOfType(store, "run_resumed");
    ASSERT_EQ(resumed.size(), 1u);
    EXPECT_EQ(resumed[0].iteration, interruptedAt);
    EXPECT_EQ(resumed[0].payload.at("iteration"), interruptedAt);
    EXPECT_EQ(resumed[0].payload.at("samples"), 100000000);
    EXPECT_EQ(eventsOfType(store, "run_started").size(), 1u);
    EXPECT_EQ(eventsOfType(store, "run_ended").size(), 2u);
    EXPECT_EQ(readJson(store / "manifest.json"), manifest);

    if (!HasFailure()) fs::remove_all(work);
}

namespace {
    // Starts a new Run in <work>/<name> from `from`.
    int startFrom(const fs::path& work, const std::string& name, const fs::path& from, std::vector<std::string> options) {
        fs::create_directories(work / name);
        options.insert(options.begin(), {"--start-from", from.string()});
        return waitFor(spawnCli(work, cliArguments(work / name, options)));
    }
}// namespace

// A new Run from a finished Run's last checkpoint: its own store, manifest and seed, iterations from 0, and a starting
// State that is the checkpoint's.
TEST(CliCheckpoint, StartFromMakesANewRunFromACheckpointsState) {
    const auto work   = makeWork("start-from");
    const auto source = work / "out";
    ASSERT_EQ(runCli(work, SHORT_RUN), 0) << consoleOutput(work);
    const auto sourceManifest = readJson(source / "manifest.json");

    ASSERT_EQ(startFrom(work, "from", source / "checkpoint.cbor", {"--seed", "9", "-b", "0", "-s", "10", "-t", "5", "--window", "10"}), 0) << consoleOutput(work);
    const auto store = work / "from";

    const auto manifest = readJson(store / "manifest.json");
    EXPECT_NE(manifest.at("run_id"), sourceManifest.at("run_id"));
    EXPECT_EQ(manifest.at("seed"), 9);
    expectEveryClosedWindowReadsBack(store);
    EXPECT_EQ(traceIterations(store), zeroTo(10, 5));
    EXPECT_EQ(eventsOfType(store, "run_started").size(), 1u);
    EXPECT_EQ(eventsOfType(store, "run_resumed").size(), 0u);
    EXPECT_EQ(scalarAt(store, "mean_coi", 0), scalarAt(source, "mean_coi", 30));
    EXPECT_NEAR(scalarAt(store, "likelihood", 0), scalarAt(source, "likelihood", 30), 1e-6);
    // The source Run is untouched.
    EXPECT_EQ(readJson(source / "manifest.json"), sourceManifest);

    if (!HasFailure()) fs::remove_all(work);
}

// A new Run's burn-in counts from 0 even when started from a checkpoint far into another Run, so its ladder adapts
// on the Sampling iterations a fresh Run's does: the 16th, then every 5th.
TEST(CliCheckpoint, AStartedRunsLadderAdaptsAsAFreshRunsDoes) {
    const auto work   = makeWork("start-from-adapts");
    const auto source = work / "out";
    auto threeChains  = SHORT_RUN;
    threeChains.insert(threeChains.end(), {"-n", "3"});
    ASSERT_EQ(runCli(work, threeChains), 0) << consoleOutput(work);

    ASSERT_EQ(startFrom(work, "from", source / "checkpoint.cbor", {"--seed", "9", "-n", "3", "-b", "22", "-s", "5", "-t", "1", "--window", "10"}), 0) << consoleOutput(work);
    std::vector<std::int64_t> adaptedAt;
    for (const auto& e : eventsOfType(work / "from", "temperatures_adapted")) {
        adaptedAt.push_back(e.iteration);
    }
    EXPECT_EQ(adaptedAt, (std::vector<std::int64_t>{16, 21}));

    if (!HasFailure()) fs::remove_all(work);
}

// A Run whose checkpoint can't be written ends with the output-failure exit code, naming the checkpoint.
TEST(CliCheckpoint, AnUnwritableCheckpointExitsWithOutputFailureNamingIt) {
    const auto work  = makeWork("unwritable-checkpoint");
    const auto store = work / "out";
    fs::create_directories(store / "checkpoint.cbor" / "occupied");

    EXPECT_EQ(runCli(work, SHORT_RUN), 3) << consoleOutput(work);
    EXPECT_NE(consoleOutput(work).find("ERROR: cannot rename into place output file " + (store / "checkpoint.cbor").string()), std::string::npos) << consoleOutput(work);

    if (!HasFailure()) fs::remove_all(work);
}
