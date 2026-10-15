//
// The Run assembly: a Run assembled from its recipe, its inputs and, for a Run resumed or started from a checkpoint,
// that checkpoint, the one way main and these tests build one. These tests check what the assembly decides: the Run
// manifest's identity and seed, where a resumed Run starts, and the two rules of a Run started from a checkpoint: its
// burn-in counts from 0, so its ladder adapts as a fresh Run's does, and it draws from its own seed.
//

#include "gtest/gtest.h"

#include "TestRuns.h"

#include <algorithm>
#include <memory>
#include <string>
#include <variant>
#include <vector>

using namespace transmission_nets::test::runs;

namespace {
    // The iterations the ladder adapted at, from the store's events.
    std::vector<std::int64_t> ladderAdaptedAt(const MemoryStoreSink& store) {
        std::vector<std::int64_t> out;
        for (const auto& e : eventsOfType(store, "temperatures_adapted")) {
            out.push_back(e.iteration);
        }
        return out;
    }

    // A finished Run in `store` with its Checkpoints in `dir`, and its latest Checkpoint.
    checkpoint::Checkpoint finishedRun(MemoryStoreSink& store, const fs::path& dir, Recipe recipe) {
        auto run = freshRun(store, recipe);
        EXPECT_TRUE(std::holds_alternative<Completed>(run.drive(neverStop)));
        return checkpoint::load(dir / checkpoint::LATEST_FILE);
    }
}// namespace

// A new Run's manifest records the provenance the caller gives, its seed and the time it was assembled at, and the
// store gets that manifest. Its Checkpoints carry the identity the manifest gives.
TEST(RunAssemblyTest, ANewRunsManifestCarriesItsProvenanceSeedAndStartTime) {
    TempDir dir;
    MemoryStoreSink store;
    const RunPlan plan{.burnIn = 2, .samples = 3};
    auto run = freshRun(store, newRecipe(dir.path, plan));

    EXPECT_EQ(run.manifest.runId, "run-1");
    EXPECT_EQ(run.manifest.seed, 3);
    EXPECT_EQ(run.manifest.startedAt, "2026-10-11T00:00:00.000Z");
    EXPECT_EQ(run.manifest.cliArguments, PROVENANCE.cliArguments);
    EXPECT_EQ(run.manifest.gitCommit, "abc123");
    EXPECT_EQ(store.manifest().at("run_id"), "run-1");
    EXPECT_EQ(store.manifest().at("started_at"), "2026-10-11T00:00:00.000Z");
    EXPECT_TRUE(std::holds_alternative<FreshRun>(run.start));
    EXPECT_EQ(run.plan.burnIn, 2);
    EXPECT_EQ(run.plan.samples, 3);

    ASSERT_TRUE(std::holds_alternative<Completed>(run.drive(neverStop)));
    EXPECT_EQ(checkpoint::load(dir.path / checkpoint::LATEST_FILE).identity, IDENTITY);
}

// A recipe with a seed of -1 gets one drawn, and the manifest records the seed drawn, not -1.
TEST(RunAssemblyTest, ASeedOfMinusOneIsDrawnAndRecorded) {
    TempDir dir;
    MemoryStoreSink store;
    auto recipe          = newRecipe(dir.path, {.burnIn = 1, .samples = 1});
    recipe.manifest.seed = -1;
    const auto run       = freshRun(store, recipe);

    EXPECT_NE(run.manifest.seed, -1);
    EXPECT_EQ(store.manifest().at("seed"), run.manifest.seed);
}

// A resumed Run keeps the manifest it was started with, whatever provenance it is assembled with now, and starts where
// its Checkpoint left it, with the Checkpoint's Parent-set ids.
TEST(RunAssemblyTest, AResumedRunKeepsItsManifestAndStartsWhereItsCheckpointLeftIt) {
    TempDir dir;
    MemoryStoreSink store;
    const RunPlan plan{.burnIn = 5, .samples = 5};
    const auto saved = finishedRun(store, dir.path, newRecipe(dir.path, plan));

    const Provenance later{.runId = "run-2", .cliArguments = {"later"}, .gitCommit = "def456"};
    auto run = assemble(resumeRecipe(store, dir.path, {.burnIn = 5, .samples = 10}), study(), saved, std::make_unique<SharedSink>(store), unread(), later, fixedClock);

    EXPECT_EQ(run.manifest.runId, "run-1");
    EXPECT_EQ(run.manifest.cliArguments, PROVENANCE.cliArguments);
    EXPECT_EQ(run.manifest.gitCommit, "abc123");
    EXPECT_EQ(run.plan.samples, 10);
    ASSERT_TRUE(std::holds_alternative<ResumedRun>(run.start));
    EXPECT_EQ(std::get<ResumedRun>(run.start).position, saved.position);
    EXPECT_EQ(std::get<ResumedRun>(run.start).parentSetIds, saved.parentSetIds);
    EXPECT_NEAR(run.ensemble->hotChain().model->getLikelihood(), scalarAt(store, "likelihood", 10), 1e-9);
}

// A Run started from a checkpoint far into another Run counts its burn-in from 0, so its ladder adapts on the Sampling
// iterations a fresh Run's does: the 16th, then every 5th.
TEST(RunAssemblyTest, AStartedRunsLadderAdaptsAsAFreshRunsDoes) {
    TempDir sourceDir;
    MemoryStoreSink source;
    auto sourceRecipe            = newRecipe(sourceDir.path, {.burnIn = 30, .samples = 5});
    sourceRecipe.manifest.chains = 3;
    const auto saved             = finishedRun(source, sourceDir.path, sourceRecipe);
    EXPECT_EQ(ladderAdaptedAt(source), (std::vector<std::int64_t>{16, 21, 26}));

    TempDir dir;
    MemoryStoreSink store;
    auto recipe            = newRecipe(dir.path, {.burnIn = 22, .samples = 5});
    recipe.manifest.chains = 3;
    recipe.manifest.seed   = 9;
    auto run               = startedRun(store, recipe, saved);
    EXPECT_EQ(run.ensemble->exchangeState().burnInIterations, 0);
    EXPECT_EQ(run.ensemble->exchangeState().barrierIterations, 0);
    EXPECT_EQ(run.ensemble->exchangeState().temperatures, saved.ensemble.exchange.temperatures);
    ASSERT_TRUE(std::holds_alternative<Completed>(run.drive(neverStop)));

    EXPECT_EQ(ladderAdaptedAt(store), (std::vector<std::int64_t>{16, 21}));
    EXPECT_EQ(eventsOfType(store, "run_started").size(), 1u);
    EXPECT_EQ(eventsOfType(store, "run_resumed").size(), 0u);
}

// A Run started from a checkpoint begins at the checkpoint's State and draws from its own seed from there: two started
// Runs with the same seed are the same row for row, and two with different seeds part ways.
TEST(RunAssemblyTest, AStartedRunBeginsAtTheCheckpointsStateAndDrawsFromItsOwnSeed) {
    TempDir sourceDir;
    MemoryStoreSink source;
    const auto saved = finishedRun(source, sourceDir.path, newRecipe(sourceDir.path, {.burnIn = 5, .samples = 10}));

    const auto started = [&](const long seed, MemoryStoreSink& store) {
        TempDir dir;
        auto recipe          = newRecipe(dir.path, {.burnIn = 0, .samples = 10});
        recipe.manifest.seed = seed;
        auto run             = startedRun(store, recipe, saved);
        EXPECT_EQ(run.manifest.seed, seed);
        ASSERT_TRUE(std::holds_alternative<Completed>(run.drive(neverStop)));
    };
    MemoryStoreSink nine;
    MemoryStoreSink nineAgain;
    MemoryStoreSink ten;
    started(9, nine);
    started(9, nineAgain);
    started(10, ten);

    for (const auto* store : {&nine, &nineAgain, &ten}) {
        EXPECT_EQ(scalarAt(*store, "mean_coi", 0), scalarAt(source, "mean_coi", 15));
        EXPECT_NEAR(scalarAt(*store, "likelihood", 0), scalarAt(source, "likelihood", 15), 1e-6);
    }
    EXPECT_EQ(nine.table("run_scalars").columns, nineAgain.table("run_scalars").columns);
    EXPECT_NE(nine.table("run_scalars").column<double>("mean_coi"), ten.table("run_scalars").column<double>("mean_coi"));
}
