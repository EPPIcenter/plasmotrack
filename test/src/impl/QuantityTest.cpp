//
// Every Quantity the State enumerates is named once: the name is the target in a sampler's id, a Trace column of the
// run_scalars or infections table, and a Checkpoint key. These tests hold the three consumers to the enumeration.
//

#include "gtest/gtest.h"

#include "core/io/MemoryStoreSink.h"
#include "impl/model/Model/Checkpoint.h"
#include "impl/model/Model/EnsembleCheckpoint.h"
#include "impl/model/Model/Ensemble.h"
#include "impl/model/Model/Model.h"
#include "impl/model/Model/RunRecorder.h"
#include "impl/model/Model/SamplerSet.h"
#include "impl/model/Model/State.h"

#include <boost/random.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

using namespace transmission_nets::impl::Model;
using transmission_nets::core::io::MemoryStoreSink;

namespace {
    using Engine = boost::random::mt19937;

    nlohmann::json exampleInput() {
        std::ifstream in(EXAMPLE_INPUT_PATH);
        return nlohmann::json::parse(in);
    }
    const std::vector<double> idp(60, 1.0 / 60);

    const std::vector<std::string> SCALARS{"mean_coi", "mean_strains_tx", "parent_set_size_prob"};
    const std::vector<std::string> PER_INFECTION{"expected_false_positives", "expected_false_negatives", "infection_duration"};

    std::vector<std::string> joined(std::vector<std::string> head, const std::vector<std::string>& tail) {
        head.insert(head.end(), tail.begin(), tail.end());
        return head;
    }
}// namespace

TEST(QuantityTest, ScalarsAndPerInfectionValuesAreNamedInEnumerationOrder) {
    const State state(studyFromJson(exampleInput(), idp, idp), std::make_shared<Engine>(1));
    EXPECT_EQ(state.scalarNames(), SCALARS);
    EXPECT_EQ(state.perInfectionNames(), PER_INFECTION);
}

TEST(QuantityTest, EveryScalarAndPerInfectionQuantityIsASamplersTarget) {
    auto rng   = std::make_shared<Engine>(1);
    auto state = std::make_shared<State>(studyFromJson(exampleInput(), idp, idp), rng);
    auto model = std::make_shared<Model>(state);
    const SamplerSet<Model> samplers(state, model, rng);
    const auto& ids = samplers.ids();

    for (const auto& name : state->scalarNames()) {
        EXPECT_TRUE(std::ranges::any_of(ids, [&](const std::string& id) { return id.ends_with("/" + name); })) << name;
    }
    for (const auto& name : state->perInfectionNames()) {
        EXPECT_TRUE(std::ranges::any_of(ids, [&](const std::string& id) { return id.find("/" + name + "/") != std::string::npos; })) << name;
    }
}

TEST(QuantityTest, TraceColumnsAreTheQuantityNamesInEnumerationOrder) {
    auto rng   = std::make_shared<Engine>(1);
    auto state = std::make_shared<State>(studyFromJson(exampleInput(), idp, idp), rng);
    auto model = std::make_shared<Model>(state);
    const SamplerSet<Model> samplers(state, model, rng);

    auto owned = std::make_unique<MemoryStoreSink>();
    auto* sink = owned.get();
    const RunManifest manifest{.runId = "run-q", .startedAt = "2026-10-14T00:00:00Z", .seed = 1, .cliArguments = {}, .engineVersion = "1.0.0", .gitCommit = "abc", .inputSha256 = {}, .chains = 1, .burnIn = 1, .samples = 2};
    RunRecorder recorder(*state, samplers.ids(), std::move(owned), {.thin = 1, .window = 1, .parentSetFloor = 1e-9}, manifest, [] { return std::chrono::system_clock::now(); });
    const EnsembleSnapshot snapshot{.hot = {*state, *model}, .slots = {{0, 1.0, model->getPrior(), model->getLikelihood()}}, .acceptance = {std::vector<transmission_nets::core::samplers::SamplerCounts>(samplers.ids().size())}, .swaps = {}};
    recorder.recordStart(snapshot);
    for (int i = 0; i < 3; ++i) {
        recorder.recordIteration(snapshot);
    }

    const auto& tables = sink->windowTables(0);
    EXPECT_EQ(tables.at("run_scalars").names(), joined(joined({"iteration", "phase", "window"}, SCALARS), {"prior", "likelihood", "observation", "transmission", "source"}));
    EXPECT_EQ(tables.at("infections").names(), joined({"iteration", "phase", "window", "infection"}, PER_INFECTION));
}

TEST(QuantityTest, CheckpointKeysAreTheQuantityNames) {
    auto rng = std::make_shared<Engine>(1);
    const Ensemble ensemble(1, 0.5, rng, false, 1, studyFromJson(exampleInput(), idp, idp));
    const auto saved = checkpoint::capture(ensemble, {.runId = "run-q", .inputHash = "sha256:x", .engineVersion = "1.0.0"}, {.iteration = 0, .phase = checkpoint::Phase::BurnIn, .lastClosedWindow = 0, .burnInRemaining = 1, .samplingRemaining = 1}, {});

    const auto& state = saved.ensemble.replicas.at(0).state;
    std::vector<std::string> scalars, perInfection;
    for (const auto& [name, _] : state.scalars) {
        scalars.push_back(name);
    }
    for (const auto& [name, values] : state.perInfection) {
        perInfection.push_back(name);
        EXPECT_EQ(values.size(), ensemble.chains.front().state->infections.size()) << name;
    }
    EXPECT_EQ(scalars, std::vector(SCALARS.begin(), SCALARS.end()));// maps sort their keys; these happen to be in order
    auto expected = PER_INFECTION;
    std::ranges::sort(expected);
    EXPECT_EQ(perInfection, expected);
}
