//
// The Run recipe: what was given, any stored Run and any checkpoint, resolved into the recipe a Run uses or one
// refusal. These tests build each from literals and check the Recipe, or the refusal's kind.
//

#include "gtest/gtest.h"

#include "impl/model/Model/RunRecipe.h"

#include <functional>
#include <string>
#include <variant>

using namespace transmission_nets::impl::Model;

namespace {
    const std::map<std::string, std::string> INPUTS{{"input", "in"}, {"symptomatic_idp", "sym"}, {"asymptomatic_idp", "asym"}};

    // A new Run's options: an output directory and its inputs, nothing else.
    GivenOptions newRun() {
        return {.outputDir = "out", .inputSha256 = INPUTS};
    }

    Recipe recipeOf(const RecipeOrRefusal& resolved) {
        if (const auto* refusal = std::get_if<Refusal>(&resolved)) {
            ADD_FAILURE() << "refused: " << refusal->message;
            return {};
        }
        return std::get<Recipe>(resolved);
    }

    // The kind of refusal, or none if `resolved` is a Recipe.
    std::optional<RefusalKind> refusalOf(const RecipeOrRefusal& resolved) {
        if (const auto* refusal = std::get_if<Refusal>(&resolved)) {
            return refusal->kind;
        }
        return std::nullopt;
    }

    // The Run a store holds: a finished Run with the defaults, run id "run-1".
    StoredRun storedRun() {
        return {.manifest = {.runId         = "run-1",
                             .startedAt     = "2026-10-11T12:00:00Z",
                             .seed          = 42,
                             .engineVersion = "1.0.0",
                             .inputSha256   = INPUTS,
                             .chains        = 1,
                             .gradient      = 0.0,
                             .burnIn        = 10,
                             .samples       = 20},
                .options  = {.thin = 5, .window = 10, .parentSetFloor = 1e-9}};
    }

    // The engine version this build stamps on a new Run.
    std::string engineVersion() {
        return std::get<Recipe>(resolve(newRun(), std::nullopt, std::nullopt)).manifest.engineVersion;
    }

    // The checkpoint of storedRun(), taken at its end.
    checkpoint::Checkpoint storedCheckpoint() {
        checkpoint::Checkpoint saved{.identity = {.runId = "run-1", .inputHash = "in", .engineVersion = engineVersion()},
                                     .position = {.iteration = 30, .phase = Phase::Sampling, .lastClosedWindow = 2},
                                     .ensemble = {.replicas = std::vector<checkpoint::Replica>(1)}};
        saved.ensemble.exchange.temperatures = {1.0};
        return saved;
    }

    // Resuming the Run in "out", with `change` made to what was given.
    RecipeOrRefusal resume(const std::function<void(GivenOptions&)>& change = [](GivenOptions&) {},
                           const StoredRun& stored = storedRun(), const checkpoint::Checkpoint& saved = storedCheckpoint()) {
        GivenOptions given{.resume = "out", .inputSha256 = INPUTS};
        change(given);
        return resolve(given, stored, saved);
    }

    // A two-chain checkpoint of another Run, "source", to start from.
    checkpoint::Checkpoint sourceCheckpoint() {
        checkpoint::Checkpoint saved{.identity = {.runId = "source", .inputHash = "in", .engineVersion = "0.9.0"},
                                     .position = {.iteration = 30, .phase = Phase::Sampling, .lastClosedWindow = 2},
                                     .ensemble = {.replicas = std::vector<checkpoint::Replica>(2)}};
        saved.ensemble.exchange.temperatures = {1.0, 0.4};
        return saved;
    }

    // Starting a new Run in "from" from sourceCheckpoint(), with `change` made to what was given.
    RecipeOrRefusal startFrom(const std::function<void(GivenOptions&)>& change = [](GivenOptions&) {}) {
        GivenOptions given{.outputDir = "from", .startFrom = "checkpoint.cbor", .inputSha256 = INPUTS};
        change(given);
        return resolve(given, std::nullopt, sourceCheckpoint());
    }
}// namespace

TEST(RunRecipe, ANewRunTakesTheDefaultsForWhatWasNotGiven) {
    const auto recipe = recipeOf(resolve(newRun(), std::nullopt, std::nullopt));

    EXPECT_EQ(recipe.mode, Mode::New);
    EXPECT_EQ(recipe.store, "out");
    EXPECT_EQ(recipe.plan.burnIn, 5000);
    EXPECT_EQ(recipe.plan.samples, 10000);
    EXPECT_EQ(recipe.recorder.thin, 10);
    EXPECT_EQ(recipe.recorder.window, 100);
    EXPECT_EQ(recipe.recorder.parentSetFloor, 1e-9);
    EXPECT_EQ(recipe.manifest.chains, 1);
    EXPECT_EQ(recipe.manifest.gradient, 0.0);
    EXPECT_EQ(recipe.manifest.seed, -1);
    EXPECT_EQ(recipe.manifest.burnIn, 5000);
    EXPECT_EQ(recipe.manifest.samples, 10000);
    EXPECT_FALSE(recipe.manifest.nullModel);
    EXPECT_EQ(recipe.manifest.inputSha256, INPUTS);
    EXPECT_EQ(recipe.cores, 1u);
}

TEST(RunRecipe, ResumeAndStartFromCantBeGivenTogether) {
    const GivenOptions given{.resume = "out", .startFrom = "checkpoint.cbor", .inputSha256 = INPUTS};

    EXPECT_EQ(refusalOf(resolve(given, storedRun(), std::nullopt)), RefusalKind::ResumeWithStartFrom);
}

TEST(RunRecipe, ResumeTakesNoOutputDirectory) {
    const GivenOptions given{.outputDir = "elsewhere", .resume = "out", .inputSha256 = INPUTS};

    EXPECT_EQ(refusalOf(resolve(given, storedRun(), std::nullopt)), RefusalKind::ResumeWithOutputDir);
}

TEST(RunRecipe, ANewRunNeedsAnOutputDirectory) {
    const GivenOptions given{.inputSha256 = INPUTS};

    EXPECT_EQ(refusalOf(resolve(given, std::nullopt, std::nullopt)), RefusalKind::NoOutputDir);
}

TEST(RunRecipe, ANewRunsOutputDirectoryMustNotHoldARun) {
    EXPECT_EQ(refusalOf(resolve(newRun(), storedRun(), std::nullopt)), RefusalKind::OutputHoldsRun);
}

TEST(RunRecipe, ResumingAStoreThatHoldsNoRunIsRefused) {
    const GivenOptions given{.resume = "out", .inputSha256 = INPUTS};

    EXPECT_EQ(refusalOf(resolve(given, std::nullopt, std::nullopt)), RefusalKind::NothingToResume);
}

TEST(RunRecipe, GivenValuesReplaceTheDefaults) {
    auto given           = newRun();
    given.burnIn         = 2;
    given.samples        = 5;
    given.thin           = 5;
    given.window         = 10;
    given.parentSetFloor = 1e-6;
    given.chains         = 4;
    given.cores          = 2;
    given.gradient       = 0.3;
    given.seed           = 42;
    given.nullModel      = true;
    const auto recipe    = recipeOf(resolve(given, std::nullopt, std::nullopt));

    EXPECT_EQ(recipe.plan.burnIn, 2);
    EXPECT_EQ(recipe.plan.samples, 5);
    EXPECT_EQ(recipe.recorder.thin, 5);
    EXPECT_EQ(recipe.recorder.window, 10);
    EXPECT_EQ(recipe.recorder.parentSetFloor, 1e-6);
    EXPECT_EQ(recipe.manifest.chains, 4);
    EXPECT_EQ(recipe.cores, 2u);
    EXPECT_EQ(recipe.manifest.gradient, 0.3);
    EXPECT_EQ(recipe.manifest.seed, 42);
    EXPECT_TRUE(recipe.manifest.nullModel);
}

TEST(RunRecipe, ThinMustBeAtLeastOne) {
    auto given = newRun();
    given.thin = 0;

    EXPECT_EQ(refusalOf(resolve(given, std::nullopt, std::nullopt)), RefusalKind::ThinOutOfRange);
}

TEST(RunRecipe, TheWindowMustBeAPositiveMultipleOfThin) {
    for (const auto window : {12, 0, -5}) {
        auto given   = newRun();
        given.thin   = 5;
        given.window = window;

        EXPECT_EQ(refusalOf(resolve(given, std::nullopt, std::nullopt)), RefusalKind::WindowOutOfRange) << "--window " << window;
    }
}

TEST(RunRecipe, TheParentSetFloorMustBeAtLeastZeroAndBelowOne) {
    for (const auto floor : {-0.1, 1.0, 1.5}) {
        auto given           = newRun();
        given.parentSetFloor = floor;

        EXPECT_EQ(refusalOf(resolve(given, std::nullopt, std::nullopt)), RefusalKind::ParentSetFloorOutOfRange) << "--parent-set-floor " << floor;
    }
    auto given           = newRun();
    given.parentSetFloor = 0;
    EXPECT_EQ(refusalOf(resolve(given, std::nullopt, std::nullopt)), std::nullopt);
}

TEST(RunRecipe, AResumedRunKeepsItsStoredRecipe) {
    auto stored               = storedRun();
    stored.manifest.nullModel = true;
    const auto recipe         = recipeOf(resume([](GivenOptions&) {}, stored));

    EXPECT_EQ(recipe.mode, Mode::Resume);
    EXPECT_EQ(recipe.store, "out");
    EXPECT_EQ(recipe.manifest.runId, "run-1");
    EXPECT_EQ(recipe.manifest.seed, 42);
    EXPECT_TRUE(recipe.manifest.nullModel);
    EXPECT_EQ(recipe.recorder.thin, 5);
    EXPECT_EQ(recipe.recorder.window, 10);
    EXPECT_EQ(recipe.plan.burnIn, 10);
    EXPECT_EQ(recipe.plan.samples, 20);
}

TEST(RunRecipe, ResumingWithTheSameOptionsGivenAgainIsAccepted) {
    const auto resolved = resume([](GivenOptions& g) {
        g.chains         = 1;
        g.thin           = 5;
        g.window         = 10;
        g.burnIn         = 10;
        g.gradient       = 0.0;
        g.parentSetFloor = 1e-9;
        g.seed           = 42;
        g.cores          = 4;
    });

    EXPECT_EQ(refusalOf(resolved), std::nullopt);
    EXPECT_EQ(recipeOf(resolved).cores, 4u);
}

TEST(RunRecipe, ResumingWithAnotherRecipeIsRefused) {
    const std::vector<std::pair<std::string, std::function<void(GivenOptions&)>>> changes{
            {"-n", [](GivenOptions& g) { g.chains = 2; }},
            {"-t", [](GivenOptions& g) { g.thin = 10; }},
            {"--window", [](GivenOptions& g) { g.window = 20; }},
            {"-b", [](GivenOptions& g) { g.burnIn = 11; }},
            {"-g", [](GivenOptions& g) { g.gradient = 0.5; }},
            {"--parent-set-floor", [](GivenOptions& g) { g.parentSetFloor = 1e-6; }},
            {"--seed", [](GivenOptions& g) { g.seed = 43; }},
            {"--null-model", [](GivenOptions& g) { g.nullModel = true; }},
    };
    for (const auto& [option, change] : changes) {
        EXPECT_EQ(refusalOf(resume(change)), RefusalKind::RecipeChanged) << option;
    }
}

TEST(RunRecipe, ResumingARunStartedFromACheckpointTakesNoGradient) {
    auto stored              = storedRun();
    stored.manifest.gradient = std::nullopt;

    EXPECT_EQ(refusalOf(resume([](GivenOptions&) {}, stored)), std::nullopt);
    EXPECT_EQ(refusalOf(resume([](GivenOptions& g) { g.gradient = 0.3; }, stored)), RefusalKind::RecipeChanged);
}

TEST(RunRecipe, ResumingMayGrowSamplesButNotShrinkThem) {
    EXPECT_EQ(recipeOf(resume([](GivenOptions& g) { g.samples = 30; })).plan.samples, 30);
    EXPECT_EQ(recipeOf(resume([](GivenOptions& g) { g.samples = 20; })).plan.samples, 20);
    EXPECT_EQ(refusalOf(resume([](GivenOptions& g) { g.samples = 19; })), RefusalKind::SamplesShrunk);
}

// A Run resumed in burn-in plans the sampling iterations its checkpoint has left; an earlier resume may have grown them.
TEST(RunRecipe, AResumedRunPlansTheSamplesItsCheckpointPlanned) {
    auto saved     = storedCheckpoint();
    saved.position = {.iteration = 4, .phase = Phase::BurnIn, .burnInRemaining = 6, .samplingRemaining = 25};

    EXPECT_EQ(recipeOf(resume([](GivenOptions&) {}, storedRun(), saved)).plan.samples, 25);
    EXPECT_EQ(refusalOf(resume([](GivenOptions& g) { g.samples = 20; }, storedRun(), saved)), RefusalKind::SamplesShrunk);

    saved.position = {.iteration = 15, .phase = Phase::Sampling, .lastClosedWindow = 1, .samplingRemaining = 25};
    EXPECT_EQ(recipeOf(resume([](GivenOptions&) {}, storedRun(), saved)).plan.samples, 30);
}

TEST(RunRecipe, ResumingNeedsTheStoresCheckpointToBeOfItsRun) {
    auto saved              = storedCheckpoint();
    saved.identity.runId    = "run-2";

    EXPECT_EQ(refusalOf(resume([](GivenOptions&) {}, storedRun(), saved)), RefusalKind::RunIdMismatch);
}

TEST(RunRecipe, ResumingNeedsTheSameIdps) {
    for (const auto* idp : {"symptomatic_idp", "asymptomatic_idp"}) {
        EXPECT_EQ(refusalOf(resume([&](GivenOptions& g) { g.inputSha256[idp] = "other"; })), RefusalKind::IdpChanged) << idp;
    }
}

TEST(RunRecipe, ResumingNeedsTheSameInputEngineVersionAndChainCount) {
    EXPECT_EQ(refusalOf(resume([](GivenOptions& g) { g.inputSha256["input"] = "other"; })), RefusalKind::InputChanged);

    auto otherVersion                   = storedCheckpoint();
    otherVersion.identity.engineVersion = "0.0.0-other";
    EXPECT_EQ(refusalOf(resume([](GivenOptions&) {}, storedRun(), otherVersion)), RefusalKind::EngineVersionChanged);

    auto twoChains     = storedCheckpoint();
    twoChains.ensemble.replicas = std::vector<checkpoint::Replica>(2);
    EXPECT_EQ(refusalOf(resume([](GivenOptions&) {}, storedRun(), twoChains)), RefusalKind::ChainCountChanged);
}

TEST(RunRecipe, StartFromMakesANewRunWithItsOwnOptions) {
    const auto recipe = recipeOf(startFrom([](GivenOptions& g) {
        g.seed   = 9;
        g.burnIn = 5;
    }));

    EXPECT_EQ(recipe.mode, Mode::StartFrom);
    EXPECT_EQ(recipe.store, "from");
    EXPECT_EQ(recipe.manifest.seed, 9);
    EXPECT_EQ(recipe.plan.burnIn, 5);
    EXPECT_EQ(recipe.plan.samples, 10000);
}

TEST(RunRecipe, StartFromTakesItsChainCountFromTheCheckpoint) {
    EXPECT_EQ(recipeOf(startFrom()).manifest.chains, 2);
    EXPECT_EQ(recipeOf(startFrom([](GivenOptions& g) { g.chains = 2; })).manifest.chains, 2);
}

// The ladder is the checkpoint's, by temperature slot, so the manifest records where it came from instead of -g.
TEST(RunRecipe, StartFromIgnoresTheGradientWithANotice) {
    const auto recipe = recipeOf(startFrom([](GivenOptions& g) { g.gradient = 0.3; }));

    EXPECT_EQ(recipe.manifest.gradient, std::nullopt);
    EXPECT_EQ(recipe.manifest.startedFrom, (nlohmann::json{{"run_id", "source"}, {"engine_version", "0.9.0"}, {"temperatures", {1.0, 0.4}}}));
    EXPECT_EQ(recipe.notices.size(), 1u);
    EXPECT_TRUE(recipeOf(startFrom()).notices.empty());
}

// A Run started from a checkpoint may test other IDPs on a burned-in State; its manifest records the ones it used.
TEST(RunRecipe, StartFromAcceptsOtherIdpsAndRecordsTheirHashes) {
    const auto recipe = recipeOf(startFrom([](GivenOptions& g) {
        g.inputSha256["symptomatic_idp"]  = "sym-2";
        g.inputSha256["asymptomatic_idp"] = "asym-2";
    }));

    EXPECT_EQ(recipe.manifest.inputSha256.at("symptomatic_idp"), "sym-2");
    EXPECT_EQ(recipe.manifest.inputSha256.at("asymptomatic_idp"), "asym-2");
}

TEST(RunRecipe, StartFromNeedsTheSameInputAndChainCountButNotEngineVersion) {
    EXPECT_EQ(refusalOf(startFrom([](GivenOptions& g) { g.inputSha256["input"] = "other"; })), RefusalKind::InputChanged);
    EXPECT_EQ(refusalOf(startFrom([](GivenOptions& g) { g.chains = 3; })), RefusalKind::ChainCountChanged);
    EXPECT_EQ(refusalOf(startFrom()), std::nullopt) << "the checkpoint's engine version is 0.9.0";
}

TEST(RunRecipe, ResumingOrStartingFromNeedsACheckpoint) {
    const GivenOptions resuming{.resume = "out", .inputSha256 = INPUTS};
    const GivenOptions starting{.outputDir = "from", .startFrom = "checkpoint.cbor", .inputSha256 = INPUTS};

    EXPECT_EQ(refusalOf(resolve(resuming, storedRun(), std::nullopt)), RefusalKind::UnusableCheckpoint);
    EXPECT_EQ(refusalOf(resolve(starting, std::nullopt, std::nullopt)), RefusalKind::UnusableCheckpoint);
}
