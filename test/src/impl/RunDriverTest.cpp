//
// The Run driver advancing a Run to its end or to a stop: when Checkpoints are taken, where a resumed Run picks up, and
// how a stop request and a forced end take effect. These tests assemble a Run from the example input, drive it with the
// in-memory store sink and a temporary store directory, and check the committed Windows, the Run events and the
// Checkpoint files.
//

#include "gtest/gtest.h"

#include "core/io/MemoryStoreSink.h"
#include "impl/model/Model/Checkpoint.h"
#include "impl/model/Model/Progress.h"
#include "impl/model/Model/RunDriver.h"
#include "impl/model/Model/RunRecorder.h"

#include <boost/random.hpp>
#include <nlohmann/json.hpp>

#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <thread>
#include <vector>

using namespace transmission_nets::impl::Model;
using transmission_nets::core::io::MemoryStoreSink;
using transmission_nets::core::io::Table;
using transmission_nets::core::io::WindowMarker;
namespace fs = std::filesystem;

namespace {
    using Engine = checkpoint::Engine;

    nlohmann::json exampleInput() {
        std::ifstream in(EXAMPLE_INPUT_PATH);
        return nlohmann::json::parse(in);
    }

    const std::vector<double> IDP(60, 1.0 / 60);

    struct TempDir {
        fs::path path = fs::temp_directory_path() / ("plasmotrack-run-driver-" + std::to_string(::getpid()) + "-" + ::testing::UnitTest::GetInstance()->current_test_info()->name());
        TempDir() {
            fs::remove_all(path);
            fs::create_directories(path);
        }
        ~TempDir() { fs::remove_all(path); }
    };

    // A store sink that writes to a MemoryStoreSink the test keeps, so one store outlives the Run recorders of a Run
    // and its resumptions.
    class SharedSink final : public transmission_nets::core::io::StoreSink {
    public:
        explicit SharedSink(MemoryStoreSink& store) : store_(store) {}
        void writeManifest(const nlohmann::json& manifest) override { store_.writeManifest(manifest); }
        void writeStructure(const std::string& table, const Table& rows) override { store_.writeStructure(table, rows); }
        void appendRows(const std::string& table, const Table& rows) override { store_.appendRows(table, rows); }
        void allowRowGroupEnd() override { store_.allowRowGroupEnd(); }
        void commitWindow(const WindowMarker& marker) override { store_.commitWindow(marker); }
        void discardWindowsAfter(const std::int64_t lastKept) override { store_.discardWindowsAfter(lastKept); }

    private:
        MemoryStoreSink& store_;
    };

    // What a Run is made from, beside its plan.
    struct Recipe {
        int chains            = 1;
        double gradient       = 0.5;
        long seed             = 3;
        unsigned int cores    = 1;
        RecorderOptions options{.thin = 1, .window = 5};
    };

    // A Run assembled from parts, as the CLI assembles one: an ensemble, its top-level engine, a Run recorder writing to
    // `store` and a Run driver checkpointing into `dir`.
    struct AssembledRun {
        std::shared_ptr<Engine> rng;
        std::unique_ptr<Ensemble> ensemble;
        std::unique_ptr<RunRecorder> recorder;
        std::unique_ptr<Progress> unread;
        std::unique_ptr<RunDriver> driver;
        RunStart start;
    };

    const checkpoint::Identity IDENTITY{.runId = "run-1", .inputHash = "sha256:input", .engineVersion = "1.0.0"};

    // The Run recorder's clock: fixed, so two stores of the same Run are equal row for row.
    std::chrono::system_clock::time_point fixedClock() {
        return std::chrono::sys_days{std::chrono::year{2026} / 10 / 11};
    }

    RunManifest manifest(const Recipe& recipe, const RunPlan& plan) {
        return {.runId = IDENTITY.runId, .seed = recipe.seed, .engineVersion = IDENTITY.engineVersion, .chains = recipe.chains, .gradient = recipe.gradient, .burnIn = plan.burnIn, .samples = plan.samples};
    }

    // The progress adapter in tests: records each call, with where the Run stood.
    class RecordingProgress final : public Progress {
    public:
        struct Call {
            std::string kind;
            RunPosition at;
            std::optional<WindowSummary> window{};
            std::optional<Heartbeat> beat{};
            std::optional<TemperatureAdaptation> adaptation{};
            std::string signal{};
            bool resumed = false;
        };
        std::vector<Call> calls;

        void started(const RunPosition& at, const bool resumed) override { calls.push_back({.kind = "started", .at = at, .resumed = resumed}); }
        void windowClosed(const RunPosition& at, const WindowSummary& window) override { calls.push_back({.kind = "windowClosed", .at = at, .window = window}); }
        void ladderAdapted(const RunPosition& at, const TemperatureAdaptation& adaptation) override { calls.push_back({.kind = "ladderAdapted", .at = at, .adaptation = adaptation}); }
        void heartbeat(const RunPosition& at, const Heartbeat& beat) override { calls.push_back({.kind = "heartbeat", .at = at, .beat = beat}); }
        void stopped(const RunPosition& at, const std::string& signal) override { calls.push_back({.kind = "stopped", .at = at, .signal = signal}); }
        void ended(const RunPosition& at) override { calls.push_back({.kind = "ended", .at = at}); }

        [[nodiscard]] std::vector<std::string> kinds() const {
            std::vector<std::string> out;
            for (const auto& c : calls) {
                out.push_back(c.kind);
            }
            return out;
        }
        [[nodiscard]] std::vector<Call> of(const std::string& kind) const {
            std::vector<Call> out;
            std::ranges::copy_if(calls, std::back_inserter(out), [&](const Call& c) { return c.kind == kind; });
            return out;
        }
    };

    // Where a Run assembled for a test reports its progress, and the clock its Run recorder and Run driver read.
    struct Reporting {
        Progress* progress = nullptr;// none: a recording nobody reads
        RunRecorder::Clock clock = fixedClock;
    };

    void assembleDriver(AssembledRun& run, MemoryStoreSink& store, const fs::path& dir, const Recipe& recipe, const RunPlan& plan, const Reporting& reporting) {
        if (reporting.progress == nullptr) {
            run.unread   = std::make_unique<RecordingProgress>();
        }
        run.recorder = std::make_unique<RunRecorder>(*run.ensemble->hotChain().state, run.ensemble->samplerIds(), std::make_unique<SharedSink>(store), recipe.options, manifest(recipe, plan), reporting.clock);
        run.driver   = std::make_unique<RunDriver>(*run.ensemble, *run.rng, *run.recorder, dir, IDENTITY, reporting.progress ? *reporting.progress : *run.unread, reporting.clock);
    }

    // A new Run.
    AssembledRun freshRun(MemoryStoreSink& store, const fs::path& dir, const Recipe& recipe, const RunPlan& plan, const Reporting& reporting = {}) {
        AssembledRun run;
        run.rng      = std::make_shared<Engine>(recipe.seed);
        run.ensemble = std::make_unique<Ensemble>(recipe.chains, recipe.gradient, run.rng, false, recipe.cores, exampleInput(), IDP, IDP);
        run.start    = FreshRun{};
        assembleDriver(run, store, dir, recipe, plan, reporting);
        return run;
    }

    // The Run whose store is `store` and `dir`, resumed from the latest Checkpoint in `dir`.
    AssembledRun resumedRun(MemoryStoreSink& store, const fs::path& dir, const Recipe& recipe, const RunPlan& plan, const Reporting& reporting = {}) {
        const auto saved = checkpoint::load(dir / "checkpoint.cbor");
        AssembledRun run;
        run.rng      = std::make_shared<Engine>();
        run.ensemble = checkpoint::restore(saved, run.rng, exampleInput(), IDP, IDP, false, recipe.cores);
        run.start    = ResumedRun{.position = saved.position, .parentSetIds = saved.parentSetIds};
        assembleDriver(run, store, dir, recipe, plan, reporting);
        return run;
    }

    // Never asks the Run to stop.
    std::optional<std::string> neverStop() {
        return std::nullopt;
    }

    // One row of the events table, its payload parsed.
    struct Event {
        std::int64_t iteration;
        std::string type;
        nlohmann::json payload;
    };

    std::vector<Event> events(const MemoryStoreSink& store) {
        const auto t = store.table("events");
        std::vector<Event> out;
        for (std::size_t r = 0; r < t.rows(); ++r) {
            out.push_back({t.column<std::int64_t>("iteration")[r], t.column<std::string>("type")[r], nlohmann::json::parse(t.column<std::string>("payload")[r])});
        }
        return out;
    }

    std::vector<Event> eventsOfType(const MemoryStoreSink& store, const std::string& type) {
        auto out = events(store);
        std::erase_if(out, [&](const Event& e) { return e.type != type; });
        return out;
    }

    std::vector<std::int64_t> windowNumbers(const MemoryStoreSink& store) {
        std::vector<std::int64_t> out;
        for (const auto& w : store.windows()) {
            out.push_back(w.window);
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

    // A run_scalars column at one iteration.
    double scalarAt(const MemoryStoreSink& store, const std::string& column, const std::int64_t iteration) {
        const auto t          = store.table("run_scalars");
        const auto iterations = t.column<std::int64_t>("iteration");
        const auto at         = std::ranges::find(iterations, iteration);
        EXPECT_NE(at, iterations.end()) << "no row at iteration " << iteration;
        return at == iterations.end() ? 0 : t.column<double>(column)[at - iterations.begin()];
    }

    // The replica in the Hot temperature slot.
    const checkpoint::Replica& hotReplica(const checkpoint::Checkpoint& saved) {
        return saved.replicas.at(saved.exchange.replicas.front());
    }
}// namespace

// Every Window's Checkpoint holds the State that Window ends at, including the short Window burn-in end closes and the
// Window the Run's end closes. The stop query, polled between Sampling iterations, sees each Checkpoint as it is.
TEST(RunDriverTest, ACheckpointIsTakenAtEachWindowsEndWithTheStateItEndsAt) {
    TempDir dir;
    MemoryStoreSink store;
    const Recipe recipe{.chains = 2};
    const RunPlan plan{.burnIn = 13, .samples = 12};
    auto run = freshRun(store, dir.path, recipe, plan);

    std::vector<checkpoint::Checkpoint> seen;
    const auto watchCheckpoints = [&]() -> std::optional<std::string> {
        if (fs::exists(dir.path / "checkpoint.cbor")) {
            auto saved = checkpoint::load(dir.path / "checkpoint.cbor");
            if (seen.empty() or seen.back().position != saved.position) {
                seen.push_back(std::move(saved));
            }
        }
        return std::nullopt;
    };
    const auto outcome = run.driver->drive(run.start, plan, watchCheckpoints);

    EXPECT_TRUE(std::holds_alternative<Completed>(outcome));
    // Burn-in: Windows 0-1 full, 2 closed short by burn-in end at 13; sampling: 3 and 4 full, 5 closed by the Run's end.
    EXPECT_EQ(windowNumbers(store), zeroTo(5));
    seen.push_back(checkpoint::load(dir.path / "checkpoint.cbor"));
    std::vector<std::int64_t> checkpointed;
    for (const auto& saved : seen) {
        SCOPED_TRACE(saved.position.iteration);
        checkpointed.push_back(saved.position.iteration);
        EXPECT_EQ(saved.identity, IDENTITY);
        EXPECT_EQ(saved.position.lastClosedWindow, static_cast<std::int64_t>(checkpointed.size()) - 1);
        EXPECT_EQ(hotReplica(saved).state.meanCoi, scalarAt(store, "mean_coi", saved.position.iteration));
    }
    EXPECT_EQ(checkpointed, (std::vector<std::int64_t>{5, 10, 13, 18, 23, 25}));
    EXPECT_EQ(seen.back().position, (RunPosition{.iteration = 25, .phase = Phase::Sampling, .lastClosedWindow = 5, .burnInRemaining = 0, .samplingRemaining = 0}));
}

// checkpoint-burn-in.cbor holds the State at burn-in end, here partway through a Window, and restores to it.
TEST(RunDriverTest, TheBurnInCheckpointHoldsTheStateAtBurnInEnd) {
    TempDir dir;
    MemoryStoreSink store;
    const Recipe recipe{.chains = 2};
    const RunPlan plan{.burnIn = 13, .samples = 12};
    auto run = freshRun(store, dir.path, recipe, plan);
    ASSERT_TRUE(std::holds_alternative<Completed>(run.driver->drive(run.start, plan, neverStop)));

    const auto burnIn = checkpoint::load(dir.path / "checkpoint-burn-in.cbor");
    EXPECT_EQ(burnIn.identity, IDENTITY);
    EXPECT_EQ(burnIn.position, (RunPosition{.iteration = 13, .phase = Phase::Sampling, .lastClosedWindow = 2, .burnInRemaining = 0, .samplingRemaining = 12}));
    EXPECT_EQ(hotReplica(burnIn).state.meanCoi, scalarAt(store, "mean_coi", 13));
    auto rng            = std::make_shared<Engine>();
    const auto restored = checkpoint::restore(burnIn, rng, exampleInput(), IDP, IDP, false, 1);
    EXPECT_NEAR(restored->hotChain().model->getLikelihood(), scalarAt(store, "likelihood", 13), 1e-9);
    EXPECT_NEAR(restored->hotChain().model->getPrior(), scalarAt(store, "prior", 13), 1e-9);
}

namespace {
    // Asks the Run to stop once it has run `iterations` Sampling iterations.
    RunDriver::StopQuery stopAfter(const std::int64_t iterations) {
        return [iterations, polls = std::int64_t{0}]() mutable -> std::optional<std::string> {
            if (polls++ == iterations) {
                return "SIGINT";
            }
            return std::nullopt;
        };
    }

    // The run_scalars Trace of `store` is that of `expected`, in every Window: the State exactly, and the log densities,
    // which a restored ensemble recomputes from the State, to rounding.
    void expectSameTrace(const MemoryStoreSink& store, const MemoryStoreSink& expected) {
        const auto actual = store.table("run_scalars");
        const auto wanted = expected.table("run_scalars");
        EXPECT_EQ(actual.column<std::int64_t>("iteration"), wanted.column<std::int64_t>("iteration"));
        EXPECT_EQ(actual.column<std::string>("phase"), wanted.column<std::string>("phase"));
        for (const auto* column : {"mean_coi", "mean_strains_tx", "parent_set_size_prob"}) {
            EXPECT_EQ(actual.column<double>(column), wanted.column<double>(column)) << column;
        }
        for (const auto* column : {"prior", "likelihood", "observation", "transmission", "source"}) {
            const auto a = actual.column<double>(column);
            const auto w = wanted.column<double>(column);
            ASSERT_EQ(a.size(), w.size()) << column;
            for (std::size_t r = 0; r < a.size(); ++r) {
                EXPECT_NEAR(a[r], w[r], 1e-9) << column << " at row " << r;
            }
        }
    }

    // Each committed Window starts where the one before it ended, from the starting State to `lastIteration`.
    void expectWindowsCover(const MemoryStoreSink& store, const std::int64_t lastIteration) {
        const auto windows = store.windows();
        ASSERT_FALSE(windows.empty());
        EXPECT_EQ(windowNumbers(store), zeroTo(windows.back().window));
        EXPECT_EQ(windows.front().firstIteration, 0);
        for (std::size_t k = 1; k < windows.size(); ++k) {
            EXPECT_EQ(windows[k].firstIteration, windows[k - 1].lastIteration + 1) << "Window " << windows[k].window;
        }
        EXPECT_EQ(windows.back().lastIteration, lastIteration);
    }

    class RunDriverStopTest : public ::testing::TestWithParam<std::int64_t> {};
}// namespace

// A stop request cuts the Run before any chosen Sampling iteration: at the start, mid Window, at a Window's end and at
// burn-in end, in either Phase. The resumed Run continues with no Window or Trace row duplicated or missing, and its
// Trace rows are the ones the Run would have written uncut.
TEST_P(RunDriverStopTest, AStoppedRunResumesWhereItStopped) {
    const std::int64_t stopAt = GetParam();
    const Recipe recipe{.chains = 2};
    const RunPlan plan{.burnIn = 13, .samples = 12};

    TempDir uncutDir;
    MemoryStoreSink uncut;
    {
        auto run = freshRun(uncut, uncutDir.path, recipe, plan);
        ASSERT_TRUE(std::holds_alternative<Completed>(run.driver->drive(run.start, plan, neverStop)));
    }

    TempDir dir;
    MemoryStoreSink store;
    {
        auto run           = freshRun(store, dir.path, recipe, plan);
        const auto outcome = run.driver->drive(run.start, plan, stopAfter(stopAt));
        ASSERT_TRUE(std::holds_alternative<Stopped>(outcome));
        EXPECT_EQ(std::get<Stopped>(outcome).signal, "SIGINT");
    }
    EXPECT_EQ(checkpoint::load(dir.path / "checkpoint.cbor").position.iteration, stopAt);
    {
        auto run = resumedRun(store, dir.path, recipe, plan);
        ASSERT_TRUE(std::holds_alternative<Completed>(run.driver->drive(run.start, plan, neverStop)));
    }

    expectWindowsCover(store, 25);
    EXPECT_EQ(store.table("run_scalars").column<std::int64_t>("iteration"), zeroTo(25));
    expectSameTrace(store, uncut);

    const auto resumed = eventsOfType(store, "run_resumed");
    ASSERT_EQ(resumed.size(), 1u);
    EXPECT_EQ(resumed[0].payload.at("iteration"), stopAt);
    EXPECT_EQ(resumed[0].payload.at("samples"), 12);
    const auto interrupted = eventsOfType(store, "interrupted");
    ASSERT_EQ(interrupted.size(), 1u);
    EXPECT_EQ(interrupted[0].iteration, stopAt);
    EXPECT_EQ(interrupted[0].payload, (nlohmann::json{{"signal", "SIGINT"}, {"second", false}}));
    EXPECT_EQ(eventsOfType(store, "run_started").size(), 1u);
    EXPECT_EQ(eventsOfType(store, "burn_in_ended").size(), 1u);
    EXPECT_EQ(eventsOfType(store, "run_ended").size(), 2u);
}

// Before the first Sampling iteration; mid Window and at a Window's end in burn-in; at burn-in end; mid Window, at a
// Window's end and before the last Sampling iteration in sampling.
INSTANTIATE_TEST_SUITE_P(EitherPhase, RunDriverStopTest, ::testing::Values(0, 3, 5, 13, 16, 23, 24));

// A Run killed after a Window's marker but before its Checkpoint is a Window ahead of the Checkpoint. Here the store is
// two Windows ahead: the Run finished, and its Checkpoint was set back to the one taken at burn-in end (after Window
// 0). Resuming discards every Window the Checkpoint doesn't cover before running on.
TEST(RunDriverTest, AResumedRunDiscardsWindowsNewerThanItsCheckpoint) {
    TempDir dir;
    MemoryStoreSink store;
    const Recipe recipe{.options = {.thin = 5, .window = 10}};
    const RunPlan plan{.burnIn = 10, .samples = 20};
    {
        auto run = freshRun(store, dir.path, recipe, plan);
        ASSERT_TRUE(std::holds_alternative<Completed>(run.driver->drive(run.start, plan, neverStop)));
    }
    ASSERT_EQ(windowNumbers(store), zeroTo(2));
    fs::copy_file(dir.path / "checkpoint-burn-in.cbor", dir.path / "checkpoint.cbor", fs::copy_options::overwrite_existing);
    {
        auto run = resumedRun(store, dir.path, recipe, plan);
        ASSERT_TRUE(std::holds_alternative<Completed>(run.driver->drive(run.start, plan, neverStop)));
    }

    expectWindowsCover(store, 30);
    EXPECT_EQ(store.table("run_scalars").column<std::int64_t>("iteration"), zeroTo(30, 5));
    const auto resumed = eventsOfType(store, "run_resumed");
    ASSERT_EQ(resumed.size(), 1u);
    EXPECT_EQ(resumed[0].payload.at("iteration"), 10);
    EXPECT_EQ(resumed[0].payload.at("window"), 0);
    // The first Run's sampling Windows and its end are gone; only the resumed Run's remain.
    EXPECT_EQ(eventsOfType(store, "run_ended").size(), 1u);
    EXPECT_EQ(eventsOfType(store, "burn_in_ended").size(), 1u);
}

// A finished Run resumed with more sampling iterations appends them, in Windows that run on from the last.
TEST(RunDriverTest, ResumingAFinishedRunWithMoreSamplingIterationsExtendsIt) {
    TempDir dir;
    MemoryStoreSink store;
    const Recipe recipe{.options = {.thin = 5, .window = 10}};
    {
        const RunPlan plan{.burnIn = 10, .samples = 20};
        auto run = freshRun(store, dir.path, recipe, plan);
        ASSERT_TRUE(std::holds_alternative<Completed>(run.driver->drive(run.start, plan, neverStop)));
    }
    {
        const RunPlan plan{.burnIn = 10, .samples = 40};
        auto run = resumedRun(store, dir.path, recipe, plan);
        ASSERT_TRUE(std::holds_alternative<Completed>(run.driver->drive(run.start, plan, neverStop)));
    }

    expectWindowsCover(store, 50);
    EXPECT_EQ(store.table("run_scalars").column<std::int64_t>("iteration"), zeroTo(50, 5));
    const auto resumed = eventsOfType(store, "run_resumed");
    ASSERT_EQ(resumed.size(), 1u);
    EXPECT_EQ(resumed[0].payload.at("iteration"), 30);
    EXPECT_EQ(resumed[0].payload.at("samples"), 40);
    EXPECT_EQ(eventsOfType(store, "burn_in_ended").size(), 1u);
    const auto last = checkpoint::load(dir.path / "checkpoint.cbor").position;
    EXPECT_EQ(last.iteration, 50);
    EXPECT_EQ(last.samplingRemaining, 0);
}

// Replica exchange over four chains, with the Hot chain changing during sampling, run on one core and on four. Only the
// OpenMP build runs the chains in parallel; without it both runs are sequential and the test shows nothing.
TEST(RunDriverTest, TheStoreIsTheSameForOneCoreAndFour) {
#ifndef _OPENMP
    GTEST_SKIP() << "built without OpenMP, so the core count changes nothing; build with OpenMP to run this test";
#endif
    const RunPlan plan{.burnIn = 100, .samples = 200};
    const auto withCores = [&](const unsigned int cores, MemoryStoreSink& store) {
        TempDir dir;
        auto run = freshRun(store, dir.path, {.chains = 4, .gradient = 0.5, .seed = 42, .cores = cores, .options = {.thin = 5, .window = 50}}, plan);
        ASSERT_TRUE(std::holds_alternative<Completed>(run.driver->drive(run.start, plan, neverStop)));
    };
    MemoryStoreSink one;
    MemoryStoreSink four;
    withCores(1, one);
    withCores(4, four);

    ASSERT_EQ(one.windows(), four.windows());
    for (std::size_t k = 0; k < one.windows().size(); ++k) {
        ASSERT_EQ(one.windowTables(k).size(), four.windowTables(k).size()) << "Window " << k;
        for (const auto& [name, table] : one.windowTables(k)) {
            ASSERT_TRUE(four.windowTables(k).contains(name)) << name;
            EXPECT_EQ(four.windowTables(k).at(name).columns, table.columns) << name << " in Window " << k;
        }
    }
    EXPECT_GT(one.table("swaps").rows(), 0u);
}

namespace {
    // Waits until `condition` holds.
    template<typename Condition>
    void waitUntil(Condition&& condition) {
        while (!condition()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
    }
}// namespace

// A second interrupt from another thread, while the Run is between Sampling iterations 8 and 9, ends the Run: both
// interrupts are recorded, the in-flight Window is closed without a Checkpoint, and nothing is recorded after, though
// the Run's own thread runs on. The forced end keeps the recorder for good, as the process ends just after it, so this
// runs in a death-test child.
namespace {
    // Drives a Run on its own thread and forces its end from this one between Sampling iterations 8 and 9, then exits
    // 0 if the store shows the forced end and nothing after it, or 1 naming what it shows instead.
    [[noreturn]] void forceEndMidRun(const fs::path& dir) {
        MemoryStoreSink store;
        const RunPlan plan{.burnIn = 13, .samples = 12};
        auto run = freshRun(store, dir, {.chains = 2}, plan);
        std::atomic<bool> reached{false};
        std::atomic<bool> forced{false};
        std::int64_t polls = 0;
        const auto pauseAtEight = [&]() -> std::optional<std::string> {
            if (polls++ == 8) {
                reached = true;
                waitUntil([&] { return forced.load(); });
            }
            return std::nullopt;
        };
        std::thread driving([&] { run.driver->drive(run.start, plan, pauseAtEight); });
        waitUntil([&] { return reached.load(); });
        run.driver->forceEnd("SIGINT", "SIGTERM");
        forced = true;
        std::this_thread::sleep_for(std::chrono::milliseconds(100));

        std::string problems;
        const auto e = events(store);
        std::vector<std::string> last;
        for (auto k = e.size() < 4 ? 0 : e.size() - 4; k < e.size(); ++k) {
            last.push_back(e[k].type);
        }
        if (last != std::vector<std::string>{"interrupted", "interrupted", "run_ended", "window_closed"}) {
            problems += "the Run's last events are not both interrupts, its end and the Window it closed\n";
        } else if (e[e.size() - 4].payload != nlohmann::json{{"signal", "SIGINT"}, {"second", false}} or e[e.size() - 3].payload != nlohmann::json{{"signal", "SIGTERM"}, {"second", true}}) {
            problems += "the interrupts recorded are not SIGINT then SIGTERM\n";
        }
        if (store.hasUncommittedTables()) {
            problems += "rows reached the store after the forced end\n";
        }
        if (store.windows().back().lastIteration != 8) {
            problems += "the last Window does not end at Sampling iteration 8\n";
        }
        if (checkpoint::load(dir / "checkpoint.cbor").position.lastClosedWindow >= store.windows().back().window) {
            problems += "the Window the forced end closed has a Checkpoint\n";
        }
        std::fputs(problems.c_str(), stderr);
        std::fflush(stderr);
        std::_Exit(problems.empty() ? 0 : 1);
    }
}// namespace

TEST(RunDriverDeathTest, AForcedEndFromAnotherThreadRecordsBothInterruptsAndNothingAfter) {
    GTEST_FLAG_SET(death_test_style, "threadsafe");
    TempDir dir;
    EXPECT_EXIT(
            forceEndMidRun(dir.path),
            ::testing::ExitedWithCode(0), "^$");
}

// Progress: the Run driver reports through its progress interface, here a recording fake.

// A Run reports its start, then each Window as it commits, with where the Run then stands, then its end.
TEST(RunDriverProgressTest, EachClosedWindowIsReportedOnceBetweenTheStartAndTheEnd) {
    TempDir dir;
    MemoryStoreSink store;
    RecordingProgress progress;
    const RunPlan plan{.burnIn = 13, .samples = 12};
    auto run = freshRun(store, dir.path, {}, plan, {.progress = &progress});
    ASSERT_TRUE(std::holds_alternative<Completed>(run.driver->drive(run.start, plan, neverStop)));

    // Windows 0-2 in burn-in, 2 closed short at 13; 3-5 in sampling.
    EXPECT_EQ(progress.kinds(), (std::vector<std::string>{"started", "windowClosed", "windowClosed", "windowClosed", "windowClosed", "windowClosed", "windowClosed", "ended"}));
    EXPECT_EQ(progress.calls.front().at, (RunPosition{.iteration = 0, .phase = Phase::BurnIn, .lastClosedWindow = -1, .burnInRemaining = 13, .samplingRemaining = 12}));
    EXPECT_FALSE(progress.calls.front().resumed);
    std::vector<std::int64_t> windows;
    std::vector<std::int64_t> iterations;
    for (const auto& c : progress.of("windowClosed")) {
        windows.push_back(c.window->window);
        iterations.push_back(c.at.iteration);
        EXPECT_EQ(c.at.lastClosedWindow, c.window->window);
    }
    EXPECT_EQ(windows, windowNumbers(store));
    EXPECT_EQ(iterations, (std::vector<std::int64_t>{5, 10, 13, 18, 23, 25}));
    EXPECT_EQ(progress.of("windowClosed")[2].window->phase, Phase::BurnIn);
    EXPECT_EQ(progress.of("windowClosed")[3].window->phase, Phase::Sampling);
    EXPECT_EQ(progress.calls.back().at, (RunPosition{.iteration = 25, .phase = Phase::Sampling, .lastClosedWindow = 5, .burnInRemaining = 0, .samplingRemaining = 0}));
}

namespace {
    // A clock that moves on by `tick` at every reading.
    RunRecorder::Clock tickingClock(const std::chrono::milliseconds tick) {
        return [now = fixedClock(), tick]() mutable { return now += tick; };
    }
}// namespace

// What a closed Window is reported to have held is what the store holds for it: its window_closed event, its swaps
// Block and the Hot chain's log densities on its last Trace row.
TEST(RunDriverProgressTest, AClosedWindowIsReportedAsTheStoreHoldsIt) {
    TempDir dir;
    MemoryStoreSink store;
    RecordingProgress progress;
    const RunPlan plan{.burnIn = 13, .samples = 12};
    auto run = freshRun(store, dir.path, {.chains = 3}, plan, {.progress = &progress, .clock = tickingClock(std::chrono::milliseconds{7})});
    ASSERT_TRUE(std::holds_alternative<Completed>(run.driver->drive(run.start, plan, neverStop)));

    const auto reported = progress.of("windowClosed");
    const auto closed   = eventsOfType(store, "window_closed");
    ASSERT_EQ(reported.size(), closed.size());
    const auto swaps = store.table("swaps");
    for (std::size_t k = 0; k < reported.size(); ++k) {
        const auto& w = *reported[k].window;
        SCOPED_TRACE(w.window);
        EXPECT_EQ(w.samplingIterations, closed[k].payload.at("sampling_iterations").get<std::int64_t>());
        EXPECT_EQ(w.seconds, closed[k].payload.at("seconds").get<double>());
        EXPECT_GT(w.seconds, 0);
        EXPECT_EQ(w.hotPrior, scalarAt(store, "prior", reported[k].at.iteration));
        EXPECT_EQ(w.hotLikelihood, scalarAt(store, "likelihood", reported[k].at.iteration));
        std::vector<SwapCounts> block;
        for (std::size_t r = 0; r < swaps.rows(); ++r) {
            if (swaps.column<std::int64_t>("window")[r] == w.window) {
                block.push_back({static_cast<std::uint64_t>(swaps.column<std::int64_t>("accepted")[r]), static_cast<std::uint64_t>(swaps.column<std::int64_t>("attempted")[r])});
            }
        }
        EXPECT_EQ(w.swaps, block);
        EXPECT_EQ(w.swaps.size(), 2u);
    }
}

namespace {
    using namespace std::chrono_literals;

    // A clock that stands still until the stop query moves it on: by `by[n]` before Sampling iteration n + 1.
    struct SteppedClock {
        std::shared_ptr<std::chrono::system_clock::time_point> now = std::make_shared<std::chrono::system_clock::time_point>(fixedClock());
        std::map<std::int64_t, std::chrono::seconds> by;

        [[nodiscard]] RunRecorder::Clock clock() const {
            return [now = now] { return *now; };
        }
        [[nodiscard]] RunDriver::StopQuery stepping() const {
            return [now = now, by = by, polls = std::int64_t{0}]() mutable -> std::optional<std::string> {
                if (const auto step = by.find(polls++); step != by.end()) {
                    *now += step->second;
                }
                return std::nullopt;
            };
        }
    };
}// namespace

// A Sampling iteration that takes over a minute, mid Window, is followed by one heartbeat, with where the Run stands.
TEST(RunDriverProgressTest, AHeartbeatFollowsTheSamplingIterationThatEndsAMinuteWithoutAReport) {
    TempDir dir;
    MemoryStoreSink store;
    RecordingProgress progress;
    SteppedClock time{.by = {{3, 61s}}};
    const RunPlan plan{.burnIn = 13, .samples = 12};
    auto run = freshRun(store, dir.path, {}, plan, {.progress = &progress, .clock = time.clock()});
    ASSERT_TRUE(std::holds_alternative<Completed>(run.driver->drive(run.start, plan, time.stepping())));

    const auto beats = progress.of("heartbeat");
    ASSERT_EQ(beats.size(), 1u);
    EXPECT_EQ(beats[0].at, (RunPosition{.iteration = 4, .phase = Phase::BurnIn, .lastClosedWindow = -1, .burnInRemaining = 9, .samplingRemaining = 12}));
    EXPECT_EQ(beats[0].beat->hotPrior, scalarAt(store, "prior", 4));
    EXPECT_EQ(beats[0].beat->hotLikelihood, scalarAt(store, "likelihood", 4));
    EXPECT_EQ(beats[0].beat->samplingIterations, 4);
    EXPECT_EQ(beats[0].beat->seconds, 61);
    EXPECT_EQ(progress.kinds()[1], "heartbeat");
}

// A closed Window is a report: the minute before a heartbeat counts from it.
TEST(RunDriverProgressTest, AClosedWindowPutsTheHeartbeatOff) {
    TempDir dir;
    MemoryStoreSink store;
    RecordingProgress progress;
    // Window 0 commits once Sampling iteration 6 is recorded, 100 s in; nothing is reported in the 50 s before 5 or after 6.
    SteppedClock time{.by = {{4, 50s}, {5, 50s}, {6, 50s}}};
    const RunPlan plan{.burnIn = 13, .samples = 12};
    auto run = freshRun(store, dir.path, {}, plan, {.progress = &progress, .clock = time.clock()});
    ASSERT_TRUE(std::holds_alternative<Completed>(run.driver->drive(run.start, plan, time.stepping())));

    EXPECT_TRUE(progress.of("heartbeat").empty());
}

// A stop request ends the Run with a report of where it stopped and for which signal, in place of its end. Resumed,
// the Run reports that it picks up there.
TEST(RunDriverProgressTest, AStoppedRunReportsWhereItStoppedAndResumesThere) {
    TempDir dir;
    MemoryStoreSink store;
    const RunPlan plan{.burnIn = 13, .samples = 12};
    const RunPosition stoppedAt{.iteration = 16, .phase = Phase::Sampling, .lastClosedWindow = 3, .burnInRemaining = 0, .samplingRemaining = 9};
    {
        RecordingProgress progress;
        auto run = freshRun(store, dir.path, {}, plan, {.progress = &progress});
        ASSERT_TRUE(std::holds_alternative<Stopped>(run.driver->drive(run.start, plan, stopAfter(16))));
        EXPECT_EQ(progress.kinds(), (std::vector<std::string>{"started", "windowClosed", "windowClosed", "windowClosed", "windowClosed", "stopped"}));
        EXPECT_EQ(progress.calls.back().at, stoppedAt);
        EXPECT_EQ(progress.calls.back().signal, "SIGINT");
    }
    RecordingProgress progress;
    auto run = resumedRun(store, dir.path, {}, plan, {.progress = &progress});
    ASSERT_TRUE(std::holds_alternative<Completed>(run.driver->drive(run.start, plan, neverStop)));
    EXPECT_EQ(progress.calls.front().kind, "started");
    EXPECT_TRUE(progress.calls.front().resumed);
    EXPECT_EQ(progress.calls.front().at, stoppedAt);
    EXPECT_EQ(progress.calls.back().kind, "ended");
}

// Each ladder adaptation is reported after the Sampling iteration it follows, as the store records it. The console
// prints no line for one, so unlike a closed Window it does not put the heartbeat off.
TEST(RunDriverProgressTest, EachLadderAdaptationIsReported) {
    TempDir dir;
    MemoryStoreSink store;
    RecordingProgress progress;
    const RunPlan plan{.burnIn = 22, .samples = 3};
    auto run = freshRun(store, dir.path, {.chains = 3}, plan, {.progress = &progress});
    ASSERT_TRUE(std::holds_alternative<Completed>(run.driver->drive(run.start, plan, neverStop)));

    const auto adapted = eventsOfType(store, "temperatures_adapted");
    const auto reported = progress.of("ladderAdapted");
    ASSERT_EQ(adapted.size(), 2u);
    ASSERT_EQ(reported.size(), adapted.size());
    for (std::size_t k = 0; k < adapted.size(); ++k) {
        EXPECT_EQ(reported[k].at.iteration, adapted[k].iteration);
        EXPECT_EQ(reported[k].at.phase, Phase::BurnIn);
        EXPECT_EQ(reported[k].adaptation->applied, adapted[k].payload.at("applied").get<bool>());
        EXPECT_EQ(reported[k].adaptation->temperatures, adapted[k].payload.at("temperatures").get<std::vector<double>>());
    }
}

// Adaptations come every 5 burn-in Sampling iterations, here every 50 s: more often than the heartbeat interval. With no
// Window closing in burn-in, a heartbeat still goes out every 60 s.
TEST(RunDriverProgressTest, AHeartbeatFiresInBurnInThoughTheLadderAdaptsMoreOften) {
    TempDir dir;
    MemoryStoreSink store;
    RecordingProgress progress;
    SteppedClock time;
    for (std::int64_t k = 0; k < 45; ++k) {
        time.by[k] = 10s;
    }
    const RunPlan plan{.burnIn = 40, .samples = 5};
    auto run = freshRun(store, dir.path, {.chains = 3, .options = {.thin = 1, .window = 100}}, plan, {.progress = &progress, .clock = time.clock()});
    ASSERT_TRUE(std::holds_alternative<Completed>(run.driver->drive(run.start, plan, time.stepping())));

    ASSERT_EQ(progress.of("ladderAdapted").size(), 5u);
    std::vector<std::int64_t> beatsInBurnIn;
    for (const auto& beat : progress.of("heartbeat")) {
        if (beat.at.phase == Phase::BurnIn) {
            beatsInBurnIn.push_back(beat.at.iteration);
        }
    }
    EXPECT_EQ(beatsInBurnIn, (std::vector<std::int64_t>{6, 12, 18, 24, 30, 36}));
}
