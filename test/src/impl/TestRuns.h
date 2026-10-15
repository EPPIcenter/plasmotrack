//
// Runs assembled for tests, through the Run assembly as main assembles them: the example input, a Run recipe with test
// defaults, the fresh, resumed and started Run writing to an in-memory store with its Checkpoints in a temporary
// directory, and readers of that store.
//

#ifndef TRANSMISSION_NETWORKS_TEST_TESTRUNS_H
#define TRANSMISSION_NETWORKS_TEST_TESTRUNS_H

#include "gtest/gtest.h"

#include "core/io/MemoryStoreSink.h"
#include "impl/model/Model/Checkpoint.h"
#include "impl/model/Model/EnsembleCheckpoint.h"
#include "impl/model/Model/Progress.h"
#include "impl/model/Model/RunAssembly.h"
#include "impl/model/Model/RunDriver.h"
#include "impl/model/Model/RunRecorder.h"

#include <nlohmann/json.hpp>

#include <unistd.h>

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace transmission_nets::test::runs {
    using namespace transmission_nets::impl::Model;
    using transmission_nets::core::io::MemoryStoreSink;
    using transmission_nets::core::io::Table;
    using transmission_nets::core::io::WindowMarker;
    namespace fs = std::filesystem;
    using Engine = checkpoint::Engine;

    inline nlohmann::json exampleInput() {
        std::ifstream in(EXAMPLE_INPUT_PATH);
        return nlohmann::json::parse(in);
    }

    inline const std::vector<double> IDP(60, 1.0 / 60);

    // The example input with flat IDPs, loaded once.
    inline const Inputs& inputs() {
        static const Inputs loaded{.input = exampleInput(), .symptomaticIdp = IDP, .asymptomaticIdp = IDP};
        return loaded;
    }

    struct TempDir {
        fs::path path = fs::temp_directory_path() / ("plasmotrack-run-" + std::to_string(::getpid()) + "-" + ::testing::UnitTest::GetInstance()->current_test_info()->name());
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

    // The identity every test Run's Checkpoints carry, and the provenance its manifest records.
    inline const checkpoint::Identity IDENTITY{.runId = "run-1", .inputHash = "sha256:input", .engineVersion = "1.0.0"};
    inline const Provenance PROVENANCE{.runId = IDENTITY.runId, .cliArguments = {"transmission_networks_model", "example_input.json", "out"}, .gitCommit = "abc123"};

    // The clock of a test Run: fixed, so two stores of the same Run are equal row for row.
    inline std::chrono::system_clock::time_point fixedClock() {
        return std::chrono::sys_days{std::chrono::year{2026} / 10 / 11};
    }

    // Progress nobody reads.
    inline Progress& unread() {
        class Unread final : public Progress {
            void started(const RunPosition&, bool) override {}
            void windowClosed(const RunPosition&, const WindowSummary&) override {}
            void ladderAdapted(const RunPosition&, const TemperatureAdaptation&) override {}
            void heartbeat(const RunPosition&, const Heartbeat&) override {}
            void stopped(const RunPosition&, const std::string&) override {}
            void ended(const RunPosition&) override {}
        };
        static Unread progress;
        return progress;
    }

    // Where a test Run reports its progress, and the clock its Run recorder and Run driver read.
    struct Reporting {
        Progress* progress       = nullptr;// none: nobody reads it
        RunRecorder::Clock clock = fixedClock;
    };

    // A new Run's recipe with test defaults: one chain at gradient 0.5 from seed 3, thin 1, Windows of 5, its
    // Checkpoints into `dir`. A test changes what it tests.
    inline Recipe newRecipe(const fs::path& dir, const RunPlan& plan) {
        Recipe recipe{.mode = Mode::New, .store = dir, .recorder = {.thin = 1, .window = 5}, .plan = plan, .cores = 1};
        recipe.manifest = {.seed          = 3,
                           .engineVersion = IDENTITY.engineVersion,
                           .inputSha256   = {{input_names::INPUT, IDENTITY.inputHash}},
                           .chains        = 1,
                           .gradient      = 0.5,
                           .burnIn        = plan.burnIn,
                           .samples       = plan.samples};
        return recipe;
    }

    // The recipe of the Run in `store` resumed with `plan`, as the Run recipe resolves a resumption: the stored Run
    // manifest and recorder options, with its Checkpoints in `dir`.
    inline Recipe resumeRecipe(const MemoryStoreSink& store, const fs::path& dir, const RunPlan& plan) {
        const auto stored = storedRunFromJson(store.manifest());
        return {.mode = Mode::Resume, .store = dir, .manifest = stored.manifest, .recorder = stored.options, .plan = plan, .cores = 1};
    }

    // A new Run.
    inline AssembledRun freshRun(MemoryStoreSink& store, const Recipe& recipe, const Reporting& reporting = {}) {
        return assemble(recipe, inputs(), std::nullopt, std::make_unique<SharedSink>(store), reporting.progress ? *reporting.progress : unread(), PROVENANCE, reporting.clock);
    }

    // The Run whose store is `store` and `dir`, resumed from the latest Checkpoint in `dir`.
    inline AssembledRun resumedRun(MemoryStoreSink& store, const fs::path& dir, const RunPlan& plan, const Reporting& reporting = {}) {
        return assemble(resumeRecipe(store, dir, plan), inputs(), checkpoint::load(dir / checkpoint::LATEST_FILE), std::make_unique<SharedSink>(store),
                        reporting.progress ? *reporting.progress : unread(), PROVENANCE, reporting.clock);
    }

    // A new Run started from `saved`, as `recipe` describes it.
    inline AssembledRun startedRun(MemoryStoreSink& store, Recipe recipe, const checkpoint::Checkpoint& saved, const Reporting& reporting = {}) {
        recipe.mode              = Mode::StartFrom;
        recipe.manifest.gradient = std::nullopt;
        return assemble(recipe, inputs(), saved, std::make_unique<SharedSink>(store), reporting.progress ? *reporting.progress : unread(), PROVENANCE, reporting.clock);
    }

    // Never asks the Run to stop.
    inline std::optional<std::string> neverStop() {
        return std::nullopt;
    }

    // One row of the events table, its payload parsed.
    struct Event {
        std::int64_t iteration;
        std::string type;
        nlohmann::json payload;
    };

    inline std::vector<Event> events(const MemoryStoreSink& store) {
        const auto t = store.table("events");
        std::vector<Event> out;
        for (std::size_t r = 0; r < t.rows(); ++r) {
            out.push_back({t.column<std::int64_t>("iteration")[r], t.column<std::string>("type")[r], nlohmann::json::parse(t.column<std::string>("payload")[r])});
        }
        return out;
    }

    inline std::vector<Event> eventsOfType(const MemoryStoreSink& store, const std::string& type) {
        auto out = events(store);
        std::erase_if(out, [&](const Event& e) { return e.type != type; });
        return out;
    }

    // A run_scalars column at one iteration.
    inline double scalarAt(const MemoryStoreSink& store, const std::string& column, const std::int64_t iteration) {
        const auto t          = store.table("run_scalars");
        const auto iterations = t.column<std::int64_t>("iteration");
        const auto at         = std::ranges::find(iterations, iteration);
        EXPECT_NE(at, iterations.end()) << "no row at iteration " << iteration;
        return at == iterations.end() ? 0 : t.column<double>(column)[at - iterations.begin()];
    }
}// namespace transmission_nets::test::runs

#endif//TRANSMISSION_NETWORKS_TEST_TESTRUNS_H
