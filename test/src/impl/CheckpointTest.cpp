//
// A checkpoint saves everything a Run carries between Sampling iterations and restores it into a new ensemble: every
// chain's State, the temperature slots, every RNG engine, every sampler's counts and adaptation state, the Run
// position, the Parent-set id dictionary and the Run's identity. These tests go through capture, save, load and restore
// only.
//

#include "gtest/gtest.h"

#include "impl/model/Model/Checkpoint.h"
#include "impl/model/Model/EnsembleCheckpoint.h"

#include <boost/random.hpp>
#include <nlohmann/json.hpp>

#include <unistd.h>

#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace transmission_nets::impl::Model;
namespace fs = std::filesystem;

namespace {
    using Engine = boost::random::mt19937;

    nlohmann::json exampleInput() {
        std::ifstream in(EXAMPLE_INPUT_PATH);
        return nlohmann::json::parse(in);
    }

    const std::vector<double> idp(60, 1.0 / 60);

    struct Run {
        std::shared_ptr<Engine> rng;
        std::unique_ptr<Ensemble> ensemble;
    };

    // An ensemble past the start of sampler adaptation (iteration 20), with swaps between its temperature slots.
    Run burnedIn(const int numChains, const unsigned int seed, const int iterations = 40) {
        auto rng      = std::make_shared<Engine>(seed);
        auto ensemble = std::make_unique<Ensemble>(numChains, 0.5, rng, false, 1, exampleInput(), idp, idp);
        for (int ii = 0; ii < iterations; ++ii) {
            ensemble->burnin();
        }
        return {rng, std::move(ensemble)};
    }

    checkpoint::Identity identity() {
        return {.runId = "run-a", .inputHash = "sha256:input", .engineVersion = "1.0.0"};
    }

    checkpoint::RunPosition position() {
        return {.iteration = 40, .phase = checkpoint::Phase::BurnIn, .lastClosedWindow = 0, .burnInRemaining = 60, .samplingRemaining = 500};
    }

    struct TempDir {
        fs::path path = fs::temp_directory_path() / ("plasmotrack-checkpoint-" + std::to_string(::getpid()) + "-" + ::testing::UnitTest::GetInstance()->current_test_info()->name());
        TempDir() { fs::create_directories(path); }
        ~TempDir() { fs::remove_all(path); }
    };

    Run restoredFrom(const fs::path& file, const unsigned int seed = 99) {
        auto rng = std::make_shared<Engine>(seed);
        auto ensemble = checkpoint::restore(checkpoint::load(file), rng, exampleInput(), idp, idp, false, 1);
        return {rng, std::move(ensemble)};
    }
}// namespace

TEST(CheckpointTest, EveryRngEngineContinuesWithTheSameDraws) {
    TempDir dir;
    const auto file = dir.path / "checkpoint.cbor";
    auto original = burnedIn(3, 11);

    checkpoint::save(checkpoint::capture(*original.ensemble, identity(), position(), {}), file);
    auto restored = restoredFrom(file);

    ASSERT_EQ(restored.ensemble->chains.size(), 3u);
    for (size_t replica = 0; replica < 3; ++replica) {
        for (int draw = 0; draw < 5; ++draw) {
            EXPECT_EQ((*restored.ensemble->chains[replica].r)(), (*original.ensemble->chains[replica].r)()) << "replica " << replica;
        }
    }
    for (int draw = 0; draw < 5; ++draw) {
        EXPECT_EQ((*restored.rng)(), (*original.rng)());
    }
}

TEST(CheckpointTest, RestoredChainsHaveTheSameStateAndLikelihood) {
    TempDir dir;
    const auto file = dir.path / "checkpoint.cbor";
    auto original = burnedIn(3, 12);
    const auto saved = checkpoint::capture(*original.ensemble, identity(), position(), {});

    checkpoint::save(saved, file);
    auto restored = restoredFrom(file);

    const auto again = checkpoint::capture(*restored.ensemble, identity(), position(), {});
    for (size_t replica = 0; replica < 3; ++replica) {
        EXPECT_EQ(again.ensemble.replicas[replica].state, saved.ensemble.replicas[replica].state) << "replica " << replica;
    }
    for (size_t replica = 0; replica < 3; ++replica) {
        auto& was = *original.ensemble->chains[replica].model;
        auto& is  = *restored.ensemble->chains[replica].model;
        EXPECT_NEAR(is.getLikelihood(), was.getLikelihood(), 1e-8 * std::abs(was.getLikelihood())) << "replica " << replica;
        EXPECT_NEAR(is.getPrior(), was.getPrior(), 1e-8 * std::abs(was.getPrior())) << "replica " << replica;
    }
}

TEST(CheckpointTest, EachTemperatureSlotKeepsItsReplicaTemperatureAndSwapBarriers) {
    TempDir dir;
    const auto file = dir.path / "checkpoint.cbor";
    auto original = burnedIn(4, 15, 60);
    for (int ii = 0; ii < 15; ++ii) {
        original.ensemble->burnin();
    }
    ASSERT_NE(original.ensemble->exchangeState().replicas, (std::vector<size_t>{0, 1, 2, 3})) << "pick a seed whose replicas have swapped";
    auto at      = position();
    at.iteration = 75;// the burn-in iterations run, which restore takes the ladder's schedule from

    checkpoint::save(checkpoint::capture(*original.ensemble, identity(), at, {}), file);
    auto restored = restoredFrom(file);

    EXPECT_EQ(restored.ensemble->exchangeState(), original.ensemble->exchangeState());
    for (size_t replica = 0; replica < 4; ++replica) {
        EXPECT_EQ(restored.ensemble->chains[replica].model->getTemperature(), original.ensemble->chains[replica].model->getTemperature());
    }
}

TEST(CheckpointTest, EverySamplerKeepsItsCountsAndAdaptationState) {
    TempDir dir;
    const auto file = dir.path / "checkpoint.cbor";
    auto original = burnedIn(2, 14, 60);

    checkpoint::save(checkpoint::capture(*original.ensemble, identity(), position(), {}), file);
    auto restored = restoredFrom(file);

    for (size_t replica = 0; replica < 2; ++replica) {
        const auto& was = *original.ensemble->chains[replica].sampler;
        const auto& is  = *restored.ensemble->chains[replica].sampler;
        ASSERT_EQ(is.ids(), was.ids());
        for (const auto& id : was.ids()) {
            EXPECT_EQ(is.sampler(id).counts(), was.sampler(id).counts()) << id;
            EXPECT_EQ(is.sampler(id).adaptationState(), was.sampler(id).adaptationState()) << id;
        }
    }
}

namespace {
    // How often each string in `saved`, keys included, appears.
    void countStrings(const nlohmann::json& saved, std::map<std::string, size_t>& counts) {
        if (saved.is_string()) {
            ++counts[saved.get_ref<const std::string&>()];
        } else if (saved.is_object()) {
            for (const auto& [key, value] : saved.items()) {
                ++counts[key];
                countStrings(value, counts);
            }
        } else if (saved.is_array()) {
            for (const auto& value : saved) {
                countStrings(value, counts);
            }
        }
    }
}// namespace

// Every chain has the same samplers, so a checkpoint names each sampler once however many chains it holds.
TEST(CheckpointTest, EachSamplerIsNamedOnceWhateverTheChainCount) {
    TempDir dir;
    const auto file = dir.path / "checkpoint.cbor";
    auto original = burnedIn(3, 21, 5);
    const auto& ids = original.ensemble->chains[0].sampler->ids();
    ASSERT_FALSE(ids.empty());

    checkpoint::save(checkpoint::capture(*original.ensemble, identity(), position(), {}), file);

    std::ifstream in(file, std::ios::binary);
    std::map<std::string, size_t> counts;
    countStrings(nlohmann::json::from_cbor(in), counts);
    for (const auto& id : ids) {
        EXPECT_EQ(counts[id], 1u) << id;
    }
}

// Every chain has the same Infections, Latent parents and Loci, so a checkpoint names each once however many chains it
// holds.
TEST(CheckpointTest, EachInfectionLatentParentAndLocusIsNamedOnceWhateverTheChainCount) {
    TempDir dir;
    const auto file = dir.path / "checkpoint.cbor";
    auto original = burnedIn(3, 22, 5);
    const auto& state = *original.ensemble->chains[0].state;
    std::vector<std::string> names;
    for (const auto& infection : state.infections) {
        names.push_back(infection->id());
    }
    for (const auto& latentParent : state.latentParents) {
        names.push_back(latentParent->id());
    }
    for (const auto& [label, locus] : state.loci) {
        names.push_back(label);
    }

    checkpoint::save(checkpoint::capture(*original.ensemble, identity(), position(), {}), file);

    std::ifstream in(file, std::ios::binary);
    std::map<std::string, size_t> counts;
    countStrings(nlohmann::json::from_cbor(in), counts);
    for (const auto& name : names) {
        EXPECT_EQ(counts[name], 1u) << name;
    }
}

TEST(CheckpointTest, IdentityRunPositionAndParentSetIdsSurvive) {
    TempDir dir;
    const auto file = dir.path / "checkpoint.cbor";
    auto original = burnedIn(1, 16, 5);
    const checkpoint::RunPosition sampling{.iteration = 1200, .phase = checkpoint::Phase::Sampling, .lastClosedWindow = 11, .burnInRemaining = 0, .samplingRemaining = 800};
    const checkpoint::ParentSetIdDictionary parentSetIds{
            {"1", {.ids = {{{}, 0}, {{"1_copy"}, 1}}, .nextId = 2}},
            {"3", {.ids = {{{"1", "2"}, 0}, {{"2"}, 4}}, .nextId = 7}},
    };

    checkpoint::save(checkpoint::capture(*original.ensemble, identity(), sampling, parentSetIds), file);
    const auto loaded = checkpoint::load(file);

    EXPECT_EQ(loaded.identity, identity());
    EXPECT_EQ(loaded.position, sampling);
    EXPECT_EQ(loaded.parentSetIds, parentSetIds);
}

TEST(CheckpointTest, SavingAgainReplacesTheCheckpointAndLeavesNoOtherFile) {
    TempDir dir;
    const auto file = dir.path / "checkpoint.cbor";
    auto original = burnedIn(1, 20, 5);
    checkpoint::RunPosition later = position();
    later.iteration = 41;

    checkpoint::save(checkpoint::capture(*original.ensemble, identity(), position(), {}), file);
    checkpoint::save(checkpoint::capture(*original.ensemble, identity(), later, {}), file);

    EXPECT_EQ(checkpoint::load(file).position, later);
    EXPECT_EQ(std::distance(fs::directory_iterator(dir.path), fs::directory_iterator()), 1);
}

TEST(CheckpointTest, LoadRefusesAnotherFormatVersion) {
    TempDir dir;
    const auto file = dir.path / "checkpoint.cbor";
    auto original = burnedIn(1, 17, 5);
    checkpoint::save(checkpoint::capture(*original.ensemble, identity(), position(), {}), file);
    std::ifstream in(file, std::ios::binary);
    auto saved = nlohmann::json::from_cbor(in);
    in.close();
    ASSERT_EQ(saved.at("format_version"), checkpoint::FORMAT_VERSION);

    saved["format_version"] = checkpoint::FORMAT_VERSION + 1;
    const auto bytes = nlohmann::json::to_cbor(saved);
    std::ofstream(file, std::ios::binary | std::ios::trunc).write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));

    EXPECT_THROW(checkpoint::load(file), Refused);
}

TEST(CheckpointTest, LoadRefusesAFileThatIsNotACheckpoint) {
    TempDir dir;
    const auto file = dir.path / "checkpoint.cbor";
    std::ofstream(file) << "not a checkpoint";

    EXPECT_THROW(checkpoint::load(file), Refused);
}

TEST(CheckpointTest, RestoreRefusesAnInputWithOtherLoci) {
    TempDir dir;
    const auto file = dir.path / "checkpoint.cbor";
    auto original = burnedIn(1, 19, 5);
    checkpoint::save(checkpoint::capture(*original.ensemble, identity(), position(), {}), file);
    std::string renamed = exampleInput().dump();
    for (size_t at = renamed.find("AS11"); at != std::string::npos; at = renamed.find("AS11", at)) {
        renamed.replace(at, 4, "AS99");
    }

    EXPECT_THROW(checkpoint::restore(checkpoint::load(file), std::make_shared<Engine>(1), nlohmann::json::parse(renamed), idp, idp, false, 1), Refused);
}

TEST(CheckpointTest, ACheckpointSavedWhileSamplingRestoresNoPendingLadderBarriers) {
    TempDir dir;
    const auto file = dir.path / "checkpoint.cbor";
    // Burn-in ends 3 iterations after the adaptation at 15, leaving barriers pending, then sampling alternates the pairs.
    auto original = burnedIn(4, 22, 18);
    for (int ii = 0; ii < 5; ++ii) {
        original.ensemble->sample();
    }
    ASSERT_GT(original.ensemble->exchangeState().barrierIterations, 0);
    const checkpoint::RunPosition sampling{.iteration = 23, .phase = checkpoint::Phase::Sampling, .lastClosedWindow = 0, .burnInRemaining = 0, .samplingRemaining = 10};
    checkpoint::save(checkpoint::capture(*original.ensemble, identity(), sampling, {}), file);

    const auto exchange = restoredFrom(file).ensemble->exchangeState();

    EXPECT_EQ(exchange.barrierIterations, 0);
    EXPECT_EQ(exchange.swapBarriers, (std::vector<double>{0.0, 0.0, 0.0}));
    EXPECT_EQ(exchange.temperatures, original.ensemble->exchangeState().temperatures);
}
