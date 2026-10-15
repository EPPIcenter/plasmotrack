//
// The Run recorder writing a store: Trace tables on the `-t` grid in both Phases, Blocks and Run events, committed per
// Window, with the Run manifest and the structure tables written at Run start. These tests drive it through its interface with the
// in-memory store sink and check the committed tables and Window markers.
//

#include "gtest/gtest.h"

#include "core/io/MemoryStoreSink.h"
#include "impl/model/Model/Model.h"
#include "impl/model/Model/State.h"
#include "impl/model/Model/RunRecorder.h"

#include <boost/random.hpp>
#include <nlohmann/json.hpp>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <unistd.h>

using namespace transmission_nets::impl::Model;
using transmission_nets::core::io::MemoryStoreSink;
using transmission_nets::core::io::Table;
using transmission_nets::core::io::WindowMarker;
using transmission_nets::core::samplers::SamplerCounts;

namespace {
    // One replica: the example input (3 Infections, 2 Loci, 3 Latent parents) with its own random start.
    struct Replica {
        std::shared_ptr<State> state;
        std::shared_ptr<Model> model;

        explicit Replica(const unsigned int seed) {
            std::ifstream in(EXAMPLE_INPUT_PATH);
            const auto input = nlohmann::json::parse(in);
            const std::vector<double> idp(60, 1.0 / 60);
            state = std::make_shared<State>(studyFromJson(input, idp, idp), std::make_shared<boost::random::mt19937>(seed));
            model = std::make_shared<Model>(state);
        }

        [[nodiscard]] HotChain hot() const { return {*state, *model}; }

        // Sets the Latent genotypes of every Infection and Latent parent, by id then Locus, as one accepted proposal.
        void setLatentGenotypes(const std::map<std::string, std::map<std::string, std::string>>& alleles) const {
            std::vector<std::pair<std::shared_ptr<transmission_nets::core::parameters::Parameter<GeneticsImpl>>, std::string>> changes;
            for (const auto& owners : {state->infections, state->latentParents}) {
                for (const auto& owner : owners) {
                    for (const auto& [locus, bits] : alleles.at(owner->id())) {
                        changes.emplace_back(owner->latentGenotype(state->loci.at(locus)), bits);
                    }
                }
            }
            for (const auto& [genotype, _] : changes) {
                genotype->saveState(1);
            }
            for (const auto& [genotype, bits] : changes) {
                genotype->setValue(GeneticsImpl(bits));
            }
            model->value();
            for (const auto& [genotype, _] : changes) {
                genotype->acceptState();
            }
        }

        // Every Infection and Latent parent carries the same Alleles, so every Parent set of every Infection is possible.
        void makeEveryParentSetPossible() const {
            const std::map<std::string, std::string> same{{"AS1", "10000000"}, {"AS11", "000000010000"}};
            setLatentGenotypes({{"1", same}, {"2", same}, {"3", same}, {"1_copy", same}, {"2_copy", same}, {"3_copy", same}});
        }
    };

    RunManifest exampleManifest() {
        return {.runId         = "run-1",
                .startedAt     = "2026-10-11T12:00:00Z",
                .seed          = 42,
                .cliArguments  = {"transmission_networks_model", "in.json", "out", "-t", "2"},
                .engineVersion = "1.0.0",
                .gitCommit     = "abc123",
                .inputSha256   = {{"input", "aa"}, {"symptomatic_idp", "bb"}, {"asymptomatic_idp", "cc"}},
                .chains        = 4,
                .gradient      = 0.5,
                .burnIn        = 5,
                .samples       = 6};
    }

    const std::vector<std::string> SAMPLER_IDS{"ConstrainedContinuousRandomWalk/mean_coi", "SALTSampler/allele_frequencies/AS1", "ParentRelinkSampler/infection/1"};

    struct RunRecorderTest : ::testing::Test {
        Replica replica{1};
        MemoryStoreSink* sink = nullptr;

        // What the caller hands over beside the Hot chain, for an ensemble of four chains. Tests change it between
        // calls to play out a Run.
        std::vector<TemperatureSlot> slots{{0, 1.0, -1.0, -10.0}, {2, 0.8, -2.0, -20.0}, {1, 0.6, -3.0, -30.0}, {3, 0.4, -4.0, -40.0}};
        // Each replica's cumulative counts, by replica then sampler id. Replica 0 is the Hot chain unless a test moves it.
        std::vector<std::vector<SamplerCounts>> acceptance = std::vector(4, std::vector<SamplerCounts>(SAMPLER_IDS.size()));
        std::vector<SwapCounts> swaps         = std::vector<SwapCounts>(3);

        // The wall clock the recorder reads; tests move it forward.
        std::chrono::system_clock::time_point now = std::chrono::sys_days{std::chrono::year{2026} / 10 / 11} + std::chrono::hours{12};

        RunRecorder makeRecorder(const std::int64_t thin, const std::int64_t window, const double parentSetFloor = 1e-9, const RunManifest& manifest = exampleManifest()) {
            auto owned = std::make_unique<MemoryStoreSink>();
            sink       = owned.get();
            return RunRecorder(*replica.state, SAMPLER_IDS, std::move(owned), {.thin = thin, .window = window, .parentSetFloor = parentSetFloor}, manifest, [this] { return now; });
        }

        [[nodiscard]] EnsembleSnapshot ensemble(const Replica& hot) const { return {.hot = hot.hot(), .slots = slots, .acceptance = acceptance, .swaps = swaps}; }
        [[nodiscard]] EnsembleSnapshot ensemble() const { return ensemble(replica); }

        // Where adaptation ends: a random walk, a SALT sampler with one variance per Allele, and a sampler that
        // doesn't adapt.
        static FinalAdaptation finalAdaptation() {
            return {.samplers     = {{.scales = {{.variance = 0.5}}, .updates = 200}, {.scales = {{.variance = 0.1}, {.variance = 0.2}, {.variance = 0.3}}, .updates = 180}, {}},
                    .temperatures = {1.0, 0.7, 0.5, 0.2}};
        }

        void wait(const double seconds) { now += std::chrono::duration_cast<std::chrono::system_clock::duration>(std::chrono::duration<double>(seconds)); }
    };

    // One row of the events table, its payload parsed.
    struct Event {
        std::int64_t iteration;
        std::string phase;
        std::int64_t window;
        std::string wallTime;
        std::string type;
        nlohmann::json payload;
    };

    std::vector<Event> events(const Table& table) {
        std::vector<Event> out;
        for (std::size_t r = 0; r < table.rows(); ++r) {
            out.push_back({table.column<std::int64_t>("iteration")[r], table.column<std::string>("phase")[r], table.column<std::int64_t>("window")[r],
                           table.column<std::string>("wall_time")[r], table.column<std::string>("type")[r], nlohmann::json::parse(table.column<std::string>("payload")[r])});
        }
        return out;
    }

    std::vector<std::string> types(const std::vector<Event>& events) {
        std::vector<std::string> out;
        for (const auto& e : events) {
            out.push_back(e.type);
        }
        return out;
    }
}// namespace

TEST_F(RunRecorderTest, ManifestHoldsTheRecipeAndTheRecordersOptions) {
    auto recorder = makeRecorder(2, 4, 1e-6);

    const auto& m = sink->manifest();
    EXPECT_EQ(m.at("run_id"), "run-1");
    EXPECT_EQ(m.at("started_at"), "2026-10-11T12:00:00Z");
    EXPECT_EQ(m.at("seed"), 42);
    EXPECT_EQ(m.at("cli_arguments"), (std::vector<std::string>{"transmission_networks_model", "in.json", "out", "-t", "2"}));
    EXPECT_EQ(m.at("engine_version"), "1.0.0");
    EXPECT_EQ(m.at("git_commit"), "abc123");
    EXPECT_EQ(m.at("input_sha256").at("symptomatic_idp"), "bb");
    EXPECT_EQ(m.at("chains"), 4);
    EXPECT_EQ(m.at("gradient"), 0.5);
    EXPECT_EQ(m.at("burn_in"), 5);
    EXPECT_EQ(m.at("samples"), 6);
    EXPECT_EQ(m.at("null_model"), false);
    EXPECT_EQ(m.at("thin"), 2);
    EXPECT_EQ(m.at("window"), 4);
    EXPECT_EQ(m.at("parent_set_floor"), 1e-6);
    EXPECT_FALSE(m.contains("started_from"));
}

// What the recorder writes into the Run manifest, it reads back, with or without a gradient.
TEST_F(RunRecorderTest, TheRunManifestRoundTripsThroughItsCodec) {
    auto manifest = exampleManifest();
    {
        auto recorder = makeRecorder(2, 4, 1e-6, manifest);
        EXPECT_EQ(storedRunFromJson(sink->manifest()), (StoredRun{.manifest = manifest, .options = {.thin = 2, .window = 4, .parentSetFloor = 1e-6}}));
    }

    manifest.gradient    = std::nullopt;
    manifest.nullModel   = true;
    manifest.startedFrom = nlohmann::json{{"run_id", "run-0"}, {"engine_version", "0.9.0"}, {"temperatures", {1.0, 0.5}}};
    auto recorder        = makeRecorder(1, 1, 0, manifest);
    EXPECT_EQ(storedRunFromJson(sink->manifest()), (StoredRun{.manifest = manifest, .options = {.thin = 1, .window = 1, .parentSetFloor = 0}}));
}

TEST_F(RunRecorderTest, ReadingAStoreFindsItsRunManifestOrNone) {
    const auto store = std::filesystem::temp_directory_path() / ("plasmotrack-stored-run-" + std::to_string(::getpid()));
    std::filesystem::create_directories(store);
    EXPECT_EQ(readStoredRun(store), std::nullopt);

    auto recorder = makeRecorder(2, 4);
    std::ofstream(store / "manifest.json") << sink->manifest().dump();
    EXPECT_EQ(readStoredRun(store), storedRunFromJson(sink->manifest()));

    std::ofstream(store / "manifest.json", std::ios::trunc) << "{\"run_id\": ";
    try {
        readStoredRun(store);
        ADD_FAILURE() << "a corrupt Run manifest was read";
    } catch (const Refused& e) {
        EXPECT_EQ(e.refusal.kind, RefusalKind::UnreadableRunManifest);
    }
    std::filesystem::remove_all(store);
}

// Allele positions run left to right in an `alleles` bitstring, the order the input lists `allele_freqs` in.
TEST_F(RunRecorderTest, StructureListsEveryAllelePositionPerLocus) {
    auto recorder = makeRecorder(1, 1);

    const auto& loci = sink->structure("loci");
    EXPECT_EQ(loci.names(), (std::vector<std::string>{"locus", "allele"}));
    ASSERT_EQ(loci.rows(), 20u);// AS1 has 8 Alleles, AS11 has 12
    EXPECT_EQ(loci.column<std::string>("locus")[0], "AS1");
    EXPECT_EQ(loci.column<std::int64_t>("allele")[7], 7);
    EXPECT_EQ(loci.column<std::string>("locus")[8], "AS11");
    EXPECT_EQ(loci.column<std::int64_t>("allele")[8], 0);
}

TEST_F(RunRecorderTest, StructureGivesEachInfectionItsLatentParentAndSymptomaticFlag) {
    auto recorder = makeRecorder(1, 1);

    const auto& infections = sink->structure("infections");
    EXPECT_EQ(infections.names(), (std::vector<std::string>{"infection", "latent_parent", "symptomatic"}));
    EXPECT_EQ(infections.column<std::string>("infection"), (std::vector<std::string>{"1", "2", "3"}));
    EXPECT_EQ(infections.column<std::string>("latent_parent"), (std::vector<std::string>{"1_copy", "2_copy", "3_copy"}));
    EXPECT_EQ(infections.column<bool>("symptomatic"), (std::vector<bool>{true, true, false}));
}

// A sampler id is `<type>/<target>`; the target may itself hold slashes.
TEST_F(RunRecorderTest, StructureSplitsEachSamplerIdIntoTypeAndTarget) {
    auto recorder = makeRecorder(1, 1);

    const auto& samplers = sink->structure("samplers");
    EXPECT_EQ(samplers.names(), (std::vector<std::string>{"sampler", "type", "target"}));
    EXPECT_EQ(samplers.column<std::string>("sampler"), SAMPLER_IDS);
    EXPECT_EQ(samplers.column<std::string>("type"), (std::vector<std::string>{"ConstrainedContinuousRandomWalk", "SALTSampler", "ParentRelinkSampler"}));
    EXPECT_EQ(samplers.column<std::string>("target"), (std::vector<std::string>{"mean_coi", "allele_frequencies/AS1", "infection/1"}));
}

TEST_F(RunRecorderTest, StartingStateIsIterationZeroOfBurnInOnEveryTrace) {
    auto recorder = makeRecorder(2, 4);
    recorder.recordStart(ensemble());
    recorder.recordEnd();

    EXPECT_EQ(sink->windows(), (std::vector<WindowMarker>{{.window = 0, .firstIteration = 0, .lastIteration = 0, .phase = "burn_in"}}));
    const std::map<std::string, std::size_t> rowsPerIteration{
            {"run_scalars", 1}, {"infections", 3}, {"ordering", 3}, {"allele_frequencies", 20}, {"latent_genotypes", 6}};
    // Latent parents are summed out of the model and have no genotypes to record.
    EXPECT_EQ(sink->table("latent_parent_genotypes").rows(), 0u);
    for (const auto& [name, rows] : rowsPerIteration) {
        const auto table = sink->table(name);
        EXPECT_EQ(table.rows(), rows) << name;
        EXPECT_EQ(table.column<std::int64_t>("iteration"), std::vector<std::int64_t>(rows, 0)) << name;
        EXPECT_EQ(table.column<std::string>("phase"), std::vector<std::string>(rows, "burn_in")) << name;
        EXPECT_EQ(table.column<std::int64_t>("window"), std::vector<std::int64_t>(rows, 0)) << name;
    }
    EXPECT_GT(sink->table("parent_sets").rows(), 0u);
}

// -t 2, --window 4, 5 burn-in and 6 sampling Sampling iterations. Burn-in is iterations 1-5, sampling 6-11.
TEST_F(RunRecorderTest, TraceRowsFallOnTheThinGridOfEachPhaseAndWindowsCloseAtSizeAndAtBurnInEnd) {
    auto recorder = makeRecorder(2, 4);
    recorder.recordStart(ensemble());
    std::vector<bool> written;
    for (int k = 0; k < 5; ++k) {
        written.push_back(recorder.recordIteration(ensemble()));
    }
    recorder.recordBurnInEnd(finalAdaptation());
    for (int k = 0; k < 6; ++k) {
        written.push_back(recorder.recordIteration(ensemble()));
    }
    recorder.recordEnd();

    EXPECT_EQ(written, (std::vector<bool>{false, true, false, true, false, false, true, false, true, false, true}));
    EXPECT_EQ(sink->windows(), (std::vector<WindowMarker>{
                                       {.window = 0, .firstIteration = 0, .lastIteration = 4, .phase = "burn_in"},
                                       {.window = 1, .firstIteration = 5, .lastIteration = 5, .phase = "burn_in"},
                                       {.window = 2, .firstIteration = 6, .lastIteration = 9, .phase = "sampling"},
                                       {.window = 3, .firstIteration = 10, .lastIteration = 11, .phase = "sampling"},
                               }));

    const auto scalars = sink->table("run_scalars");
    EXPECT_EQ(scalars.column<std::int64_t>("iteration"), (std::vector<std::int64_t>{0, 2, 4, 7, 9, 11}));
    EXPECT_EQ(scalars.column<std::string>("phase"), (std::vector<std::string>{"burn_in", "burn_in", "burn_in", "sampling", "sampling", "sampling"}));
    EXPECT_EQ(scalars.column<std::int64_t>("window"), (std::vector<std::int64_t>{0, 0, 0, 2, 2, 3}));
    // Window 1 holds no Trace row, so it hands over its Blocks and Run events only.
    std::vector<std::string> window1;
    for (const auto& [name, _] : sink->windowTables(1)) {
        window1.push_back(name);
    }
    EXPECT_EQ(window1, (std::vector<std::string>{"acceptance", "events", "swaps"}));
}

// A first interrupt lets the Sampling iteration finish; the Run then ends, which closes the in-flight Window short.
TEST_F(RunRecorderTest, FirstInterruptEndsTheRunInTheInFlightWindow) {
    auto recorder = makeRecorder(1, 10);
    recorder.recordStart(ensemble());
    recorder.recordBurnInEnd(finalAdaptation());
    for (int k = 0; k < 3; ++k) {
        wait(1.0);
        recorder.recordIteration(ensemble());
    }
    recorder.recordInterrupt("SIGINT", Interrupt::First);
    recorder.recordEnd();

    EXPECT_EQ(sink->windows(), (std::vector<WindowMarker>{
                                       {.window = 0, .firstIteration = 0, .lastIteration = 0, .phase = "burn_in"},
                                       {.window = 1, .firstIteration = 1, .lastIteration = 3, .phase = "sampling"},
                               }));
    EXPECT_EQ(sink->table("run_scalars").column<std::int64_t>("iteration"), (std::vector<std::int64_t>{0, 1, 2, 3}));
    EXPECT_FALSE(sink->hasUncommittedTables());

    const auto e = events(sink->windowTables(1).at("events"));
    ASSERT_EQ(types(e), (std::vector<std::string>{"interrupted", "run_ended", "window_closed"}));
    EXPECT_EQ(e[0].iteration, 3);
    EXPECT_EQ(e[0].payload, (nlohmann::json{{"signal", "SIGINT"}, {"second", false}}));
    EXPECT_EQ(e[1].payload, (nlohmann::json{{"reason", "interrupted"}, {"seconds", 3.0}}));
}

// A second interrupt may come while a Sampling iteration is in flight, so the Run won't reach its end: the recorder
// ends it there.
TEST_F(RunRecorderTest, SecondInterruptClosesTheInFlightWindowAndEndsTheRun) {
    auto recorder = makeRecorder(1, 10);
    recorder.recordStart(ensemble());
    recorder.recordIteration(ensemble());
    recorder.recordInterrupt("SIGINT", Interrupt::First);
    recorder.recordInterrupt("SIGTERM", Interrupt::Second);

    EXPECT_EQ(sink->windows(), (std::vector<WindowMarker>{{.window = 0, .firstIteration = 0, .lastIteration = 1, .phase = "burn_in"}}));
    EXPECT_FALSE(sink->hasUncommittedTables());
    const auto e = events(sink->table("events"));
    ASSERT_EQ(types(e), (std::vector<std::string>{"run_started", "interrupted", "interrupted", "run_ended", "window_closed"}));
    EXPECT_EQ(e[2].payload, (nlohmann::json{{"signal", "SIGTERM"}, {"second", true}}));
    EXPECT_EQ(e[3].payload.at("reason"), "interrupted");

    // The Run has ended; nothing more is recorded.
    recorder.recordEnd();
    EXPECT_EQ(sink->windows().size(), 1u);
}

namespace {
    using Members = std::set<std::string>;

    // Each (Infection, set id) → its members, failing the test if a set is listed twice.
    std::map<std::pair<std::string, std::int64_t>, Members> membersById(const Table& members) {
        std::map<std::pair<std::string, std::int64_t>, Members> out;
        const auto& infection = members.column<std::string>("infection");
        const auto& setId     = members.column<std::int64_t>("set_id");
        const auto& member    = members.column<std::string>("member");
        for (std::size_t r = 0; r < members.rows(); ++r) {
            const auto [_, inserted] = out[{infection[r], setId[r]}].insert(member[r]);
            EXPECT_TRUE(inserted) << infection[r] << " set " << setId[r] << " lists " << member[r] << " twice";
        }
        return out;
    }
}// namespace

// With no floor, Infection 3 (Candidate parents 1 and 2) has seven Parent sets, each with the Latent parent or not.
// None is the Latent parent's source term on its own: {3_copy} is the Parent set of the Latent parent alone.
TEST_F(RunRecorderTest, ParentSetsWithNoFloorCoverEveryParentSetOnce) {
    replica.makeEveryParentSetPossible();
    auto recorder = makeRecorder(1, 10, 0.0);
    recorder.recordStart(ensemble());
    recorder.recordEnd();

    const auto sets    = sink->table("parent_sets");
    const auto members = membersById(sink->table("parent_set_members"));
    std::set<Members> infection3;
    std::map<std::string, double> total;
    for (std::size_t r = 0; r < sets.rows(); ++r) {
        const auto& infection = sets.column<std::string>("infection")[r];
        total[infection] += sets.column<double>("prob")[r];
        if (infection == "3") {
            infection3.insert(members.at({infection, sets.column<std::int64_t>("set_id")[r]}));
        }
    }
    EXPECT_EQ(infection3, (std::set<Members>{{"3_copy"}, {"1"}, {"1", "3_copy"}, {"2"}, {"2", "3_copy"}, {"1", "2"}, {"1", "2", "3_copy"}}));
    for (const auto& [infection, prob] : total) {
        EXPECT_NEAR(prob, 1.0, 1e-9) << infection;
    }
    // Infection 1 has no Candidate parents, so its Latent parent alone has all the probability.
    EXPECT_EQ(members.at({"1", 0}), (Members{"1_copy"}));
}

TEST_F(RunRecorderTest, ParentSetsBelowTheFloorAreDroppedWithoutRenormalising) {
    auto unfloored = makeRecorder(1, 10, 0.0);
    unfloored.recordStart(ensemble());
    unfloored.recordEnd();
    const auto all = sink->table("parent_sets").column<double>("prob");

    auto floored = makeRecorder(1, 10, 0.01);
    floored.recordStart(ensemble());
    floored.recordEnd();
    const auto kept = sink->table("parent_sets").column<double>("prob");

    std::vector<double> expected;
    std::ranges::copy_if(all, std::back_inserter(expected), [](const double p) { return p >= 0.01; });
    ASSERT_LT(expected.size(), all.size()) << "the example needs a Parent set below the floor";
    EXPECT_EQ(kept, expected);
}

// Infection 3 and its Latent parent carry an Allele at AS1 that neither observed Infection has, so only the Latent parent
// can supply it: every Parent set without the Latent parent has probability 0. A floor of 0 keeps every non-zero row
// and no zero one.
TEST_F(RunRecorderTest, ParentSetsWithProbabilityZeroAreDroppedEvenWithNoFloor) {
    const std::map<std::string, std::string> common{{"AS1", "10000000"}, {"AS11", "000000010000"}};
    const std::map<std::string, std::string> rare{{"AS1", "00000001"}, {"AS11", "000000010000"}};
    replica.setLatentGenotypes({{"1", common}, {"2", common}, {"3", rare}, {"1_copy", common}, {"2_copy", common}, {"3_copy", rare}});

    auto recorder = makeRecorder(1, 10, 0.0);
    recorder.recordStart(ensemble());
    recorder.recordEnd();

    const auto sets    = sink->table("parent_sets");
    const auto members = membersById(sink->table("parent_set_members"));
    std::set<Members> infection3;
    for (std::size_t r = 0; r < sets.rows(); ++r) {
        const auto& infection = sets.column<std::string>("infection")[r];
        EXPECT_GT(sets.column<double>("prob")[r], 0.0) << infection << " set " << sets.column<std::int64_t>("set_id")[r];
        if (infection == "3") {
            infection3.insert(members.at({infection, sets.column<std::int64_t>("set_id")[r]}));
        }
    }
    EXPECT_EQ(infection3, (std::set<Members>{{"3_copy"}, {"1", "3_copy"}, {"2", "3_copy"}, {"1", "2", "3_copy"}}));
}

// The Hot chain alternates between two replicas, and only one of them makes every Parent set of Infection 3 possible,
// so the Parent sets written change between iterations.
TEST_F(RunRecorderTest, SetIdsAreStableAndEachNewSetIsListedInOneWindowOnly) {
    Replica other{7};
    other.makeEveryParentSetPossible();
    auto recorder = makeRecorder(1, 2, 0.0);
    recorder.recordStart(ensemble());
    for (int k = 0; k < 6; ++k) {
        recorder.recordIteration(ensemble(k % 2 == 0 ? other : replica));
    }
    recorder.recordEnd();

    const auto members = membersById(sink->table("parent_set_members"));
    std::size_t listedInWindows = 0;
    for (std::size_t w = 0; w < sink->windows().size(); ++w) {
        if (const auto found = sink->windowTables(w).find("parent_set_members"); found != sink->windowTables(w).end()) {
            listedInWindows += membersById(found->second).size();
        }
    }
    EXPECT_EQ(listedInWindows, members.size());

    // Every row's set id is known, and within one iteration an Infection never repeats a set id.
    const auto sets = sink->table("parent_sets");
    std::set<std::tuple<std::int64_t, std::string, std::int64_t>> seen;
    for (std::size_t r = 0; r < sets.rows(); ++r) {
        const auto& infection = sets.column<std::string>("infection")[r];
        const auto id         = sets.column<std::int64_t>("set_id")[r];
        EXPECT_TRUE(members.contains({infection, id})) << infection << " set " << id;
        EXPECT_TRUE(seen.insert({sets.column<std::int64_t>("iteration")[r], infection, id}).second);
    }
    // Ids are numbered from 0 per Infection in order of first appearance: Infection 3's seven sets are 0-6.
    for (std::int64_t id = 0; id < 7; ++id) {
        EXPECT_TRUE(members.contains({"3", id})) << id;
    }
    EXPECT_FALSE(members.contains({"3", 7}));
}

// #56: an open Window hands each Trace row to the sink as it is written, so rows don't wait in the recorder until the
// Window closes. A row group may end only after a whole Trace row, so each holds whole Sampling iterations in the order
// #55 set; the sink decides which of those ends it takes.
TEST_F(RunRecorderTest, AnOpenWindowHandsOverEachTraceRowAndLetsARowGroupEndOnlyAfterIt) {
    auto recorder = makeRecorder(2, 100, 0.0);
    recorder.recordStart(ensemble());
    for (int k = 0; k < 6; ++k) {
        recorder.recordIteration(ensemble());
    }
    EXPECT_TRUE(sink->windows().empty());
    EXPECT_EQ(sink->uncommittedRows("run_scalars"), 4u);
    recorder.recordIteration(ensemble());
    recorder.recordEnd();

    ASSERT_EQ(sink->windows().size(), 1u);
    EXPECT_EQ(sink->rowGroups(0, "run_scalars"), (std::vector<std::size_t>{1, 1, 1, 1}));
    const auto sets       = sink->table("parent_sets");
    const auto& iteration = sets.column<std::int64_t>("iteration");
    std::vector<std::int64_t> groupIterations;
    std::size_t first = 0;
    for (const auto rows : sink->rowGroups(0, "parent_sets")) {
        ASSERT_GT(rows, 0u);
        EXPECT_EQ(iteration.at(first), iteration.at(first + rows - 1));
        groupIterations.push_back(iteration.at(first));
        first += rows;
    }
    EXPECT_EQ(groupIterations, (std::vector<std::int64_t>{0, 2, 4, 6}));
    EXPECT_EQ(sink->rowGroups(0, "acceptance").size(), 1u);
}

// #55: an Infection's rows in one iteration come in ascending set id, not in the order the Parent sets were scored.
// Candidate parents are scored in the Infection ordering, which the sampler changes, while ids follow first
// appearance. Sorted rows repeat from one iteration to the next, so `set_id` compresses.
TEST_F(RunRecorderTest, ParentSetRowsOfAnInfectionInOneIterationAreInAscendingSetId) {
    const auto setDuration = [&](const std::size_t infection, const double duration) {
        const auto parameter = replica.state->infections[infection]->infectionDuration();
        parameter->saveState(1);
        parameter->setValue(duration);
        parameter->acceptState();
        replica.model->value();
    };
    // Infection times are observation time (234, 306, 350) minus duration (10 to start). With Infection 3 between 2
    // and 1, its only candidate parent is 2, so its sets with 2 take the first ids.
    setDuration(1, 90);
    setDuration(2, 130);
    auto recorder = makeRecorder(1, 2, 0.0);
    recorder.recordStart(ensemble());
    // Back in place, 3 gains candidate 1, which it scores before 2 but whose sets take the next ids.
    setDuration(2, 10);
    setDuration(1, 10);
    recorder.recordIteration(ensemble());
    recorder.recordEnd();

    const auto sets       = sink->table("parent_sets");
    const auto& iteration = sets.column<std::int64_t>("iteration");
    const auto& infection = sets.column<std::string>("infection");
    const auto& setId     = sets.column<std::int64_t>("set_id");
    ASSERT_GT(sets.rows(), 1u);
    for (std::size_t r = 1; r < sets.rows(); ++r) {
        if (iteration[r] == iteration[r - 1] and infection[r] == infection[r - 1]) {
            EXPECT_LT(setId[r - 1], setId[r]) << "iteration " << iteration[r] << ", Infection " << infection[r];
        }
    }
}

TEST_F(RunRecorderTest, RowsReadWhicheverChainIsHotAtTheCall) {
    Replica other{2};
    replica.state->meanCOI->initializeValue(3.5);
    other.state->meanCOI->initializeValue(7.25);
    const auto& as1 = replica.state->loci.at("AS1");
    replica.state->infections[0]->latentGenotype(as1)->initializeValue(GeneticsImpl("00010001"));
    other.state->infections[0]->latentGenotype(other.state->loci.at("AS1"))->initializeValue(GeneticsImpl("10000000"));
    auto recorder = makeRecorder(1, 10);

    recorder.recordStart(ensemble());
    recorder.recordIteration(ensemble(other));
    recorder.recordIteration(ensemble());
    recorder.recordEnd();

    EXPECT_EQ(sink->table("run_scalars").column<double>("mean_coi"), (std::vector<double>{3.5, 7.25, 3.5}));
    const auto genotypes = sink->table("latent_genotypes");
    std::vector<std::string> infection1AS1;
    for (std::size_t r = 0; r < genotypes.rows(); ++r) {
        if (genotypes.column<std::string>("infection")[r] == "1" and genotypes.column<std::string>("locus")[r] == "AS1") {
            infection1AS1.push_back(genotypes.column<std::string>("alleles")[r]);
        }
    }
    EXPECT_EQ(infection1AS1, (std::vector<std::string>{"00010001", "10000000", "00010001"}));
}

// The likelihood is the Observation process plus the Transmission process, and Parent-set size gets no column of its
// own.
TEST_F(RunRecorderTest, ProcessTermsAccountForTheLikelihood) {
    auto recorder = makeRecorder(1, 10);
    recorder.recordStart(ensemble());
    recorder.recordEnd();

    const auto scalars       = sink->table("run_scalars");
    EXPECT_EQ(scalars.names(), (std::vector<std::string>{"iteration", "phase", "window", "mean_coi", "mean_strains_tx", "parent_set_size_prob", "prior", "likelihood", "observation", "transmission", "source"}));
    const double likelihood   = scalars.column<double>("likelihood")[0];
    const double observation  = scalars.column<double>("observation")[0];
    const double transmission = scalars.column<double>("transmission")[0];
    EXPECT_TRUE(std::isfinite(likelihood));
    EXPECT_NEAR(observation + transmission, likelihood, 1e-9 * std::abs(likelihood));
    EXPECT_NE(observation, 0.0);
    EXPECT_TRUE(std::isfinite(scalars.column<double>("prior")[0]));
    // Latent parents are summed out of the Transmission process, so there is no separate Source term; the column is
    // kept for the table's layout.
    EXPECT_EQ(scalars.column<double>("source")[0], 0.0);
}

// Each call hands over every temperature slot; slot 0 is the Hot chain. Rows follow the same `-t` grid as the other
// Traces, starting state included.
TEST_F(RunRecorderTest, ChainsHoldOneRowPerTemperatureSlotOnTheThinGrid) {
    auto recorder = makeRecorder(2, 4);
    recorder.recordStart(ensemble());
    recorder.recordIteration(ensemble());
    slots[0] = {2, 1.0, -2.5, -25.0};
    slots[1] = {0, 0.8, -1.5, -15.0};
    recorder.recordIteration(ensemble());
    recorder.recordEnd();

    const auto chains = sink->table("chains");
    EXPECT_EQ(chains.names(), (std::vector<std::string>{"iteration", "phase", "window", "temperature_slot", "replica", "temperature", "prior", "likelihood"}));
    EXPECT_EQ(chains.column<std::int64_t>("iteration"), (std::vector<std::int64_t>{0, 0, 0, 0, 2, 2, 2, 2}));
    EXPECT_EQ(chains.column<std::int64_t>("temperature_slot"), (std::vector<std::int64_t>{0, 1, 2, 3, 0, 1, 2, 3}));
    EXPECT_EQ(chains.column<std::int64_t>("replica"), (std::vector<std::int64_t>{0, 2, 1, 3, 2, 0, 1, 3}));
    EXPECT_EQ(chains.column<double>("temperature"), (std::vector<double>{1.0, 0.8, 0.6, 0.4, 1.0, 0.8, 0.6, 0.4}));
    EXPECT_EQ(chains.column<double>("prior"), (std::vector<double>{-1.0, -2.0, -3.0, -4.0, -2.5, -1.5, -3.0, -4.0}));
    EXPECT_EQ(chains.column<double>("likelihood"), (std::vector<double>{-10.0, -20.0, -30.0, -40.0, -25.0, -15.0, -30.0, -40.0}));
}

// -t 2, --window 4, 5 burn-in and 6 sampling Sampling iterations, as above: Windows hold iterations 1-4, 5, 6-9 and
// 10-11. At iteration k the first sampler proposes once and accepts, the second proposes twice and accepts when k is
// even, the third proposes three times and never accepts. Slots 0 and 1 attempt a swap when k is odd and accept one
// when k is a multiple of 3. Counts start where an earlier part of the Run left them.
TEST_F(RunRecorderTest, BlocksCountEverySamplingIterationOfTheirWindow) {
    acceptance[0] = {{.accepted = 10, .proposed = 20}, {.accepted = 0, .proposed = 0}, {.accepted = 5, .proposed = 50}};
    swaps[0]      = {.accepted = 7, .attempted = 9};
    auto recorder = makeRecorder(2, 4);
    recorder.recordStart(ensemble());
    auto iterate = [&](const std::uint64_t k) {
        acceptance[0][0].accepted += 1;
        acceptance[0][0].proposed += 1;
        acceptance[0][1].accepted += k % 2 == 0;
        acceptance[0][1].proposed += 2;
        acceptance[0][2].proposed += 3;
        swaps[0].attempted += k % 2 == 1;
        swaps[0].accepted += k % 3 == 0;
        recorder.recordIteration(ensemble());
    };
    for (std::uint64_t k = 1; k <= 5; ++k) {
        iterate(k);
    }
    recorder.recordBurnInEnd(finalAdaptation());
    for (std::uint64_t k = 6; k <= 11; ++k) {
        iterate(k);
    }
    recorder.recordEnd();

    const auto blocks = sink->table("acceptance");
    EXPECT_EQ(blocks.names(), (std::vector<std::string>{"window", "phase", "sampler", "accepted", "proposed"}));
    EXPECT_EQ(blocks.column<std::int64_t>("window"), (std::vector<std::int64_t>{0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3}));
    EXPECT_EQ(blocks.column<std::string>("phase"), (std::vector<std::string>{"burn_in", "burn_in", "burn_in", "burn_in", "burn_in", "burn_in", "sampling", "sampling", "sampling", "sampling", "sampling", "sampling"}));
    std::vector<std::string> samplers;
    for (int w = 0; w < 4; ++w) {
        samplers.insert(samplers.end(), SAMPLER_IDS.begin(), SAMPLER_IDS.end());
    }
    EXPECT_EQ(blocks.column<std::string>("sampler"), samplers);
    EXPECT_EQ(blocks.column<std::int64_t>("accepted"), (std::vector<std::int64_t>{4, 2, 0, 1, 0, 0, 4, 2, 0, 2, 1, 0}));
    EXPECT_EQ(blocks.column<std::int64_t>("proposed"), (std::vector<std::int64_t>{4, 8, 12, 1, 2, 3, 4, 8, 12, 2, 4, 6}));

    const auto swapBlocks = sink->table("swaps");
    EXPECT_EQ(swapBlocks.names(), (std::vector<std::string>{"window", "phase", "temperature_pair", "accepted", "attempted"}));
    EXPECT_EQ(swapBlocks.column<std::int64_t>("window"), (std::vector<std::int64_t>{0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3}));
    EXPECT_EQ(swapBlocks.column<std::int64_t>("temperature_pair"), (std::vector<std::int64_t>{0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2}));
    EXPECT_EQ(swapBlocks.column<std::int64_t>("accepted"), (std::vector<std::int64_t>{1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0}));
    EXPECT_EQ(swapBlocks.column<std::int64_t>("attempted"), (std::vector<std::int64_t>{2, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0}));
}

// -t 1, --window 4. Every replica's first sampler makes proposals at every Sampling iteration: replica 0 one, replica 2
// ten, replica 1 a hundred, replica 3 a thousand, each accepted. Replicas 0 and 2 swap temperature slots at the end of
// iteration 2, so replica 0 is the Hot chain for iterations 1-2 and replica 2 for iterations 3-4. A swap is the last
// thing in a Sampling iteration, so the call for iteration 2 already shows replica 2 in slot 0.
TEST_F(RunRecorderTest, AcceptanceCountsOnlyTheReplicaInTheHotSlotAtEachSamplingIteration) {
    auto recorder = makeRecorder(1, 4);
    recorder.recordStart(ensemble());
    const std::vector<std::uint64_t> proposals{1, 100, 10, 1000};// by replica
    for (int k = 1; k <= 4; ++k) {
        for (std::size_t replica = 0; replica < 4; ++replica) {
            acceptance[replica][0].accepted += proposals[replica];
            acceptance[replica][0].proposed += proposals[replica];
        }
        if (k == 2) {
            std::swap(slots[0].replica, slots[1].replica);
        }
        recorder.recordIteration(ensemble());
    }
    recorder.recordEnd();

    const auto blocks = sink->table("acceptance");
    EXPECT_EQ(blocks.column<std::int64_t>("window"), (std::vector<std::int64_t>{0, 0, 0}));
    EXPECT_EQ(blocks.column<std::int64_t>("accepted"), (std::vector<std::int64_t>{22, 0, 0}));
    EXPECT_EQ(blocks.column<std::int64_t>("proposed"), (std::vector<std::int64_t>{22, 0, 0}));
}

// -t 2, --window 4, 5 burn-in and 6 sampling Sampling iterations, each taking half a second.
TEST_F(RunRecorderTest, ACompletedRunRecordsItsStartEachWindowAndItsEnd) {
    auto recorder = makeRecorder(2, 4);
    recorder.recordStart(ensemble());
    for (int k = 1; k <= 5; ++k) {
        wait(0.5);
        recorder.recordIteration(ensemble());
    }
    recorder.recordBurnInEnd(finalAdaptation());
    for (int k = 6; k <= 11; ++k) {
        wait(0.5);
        recorder.recordIteration(ensemble());
    }
    recorder.recordEnd();

    const auto table = sink->table("events");
    EXPECT_EQ(table.names(), (std::vector<std::string>{"iteration", "phase", "window", "wall_time", "type", "payload"}));
    const auto e = events(table);
    ASSERT_EQ(types(e), (std::vector<std::string>{"run_started", "window_closed", "burn_in_ended", "window_closed", "window_closed", "run_ended", "window_closed"}));

    EXPECT_EQ(e[0].iteration, 0);
    EXPECT_EQ(e[0].phase, "burn_in");
    EXPECT_EQ(e[0].window, 0);
    EXPECT_EQ(e[0].wallTime, "2026-10-11T12:00:00.000Z");
    EXPECT_EQ(e[0].payload, (nlohmann::json{{"started_at", "2026-10-11T12:00:00.000Z"}, {"burn_in", 5}, {"samples", 6}}));

    // Iterations 1-4, 5, 6-9 and 10-11; the starting state takes no time and is not a Sampling iteration. A Window
    // that fills closes when the next Sampling iteration is recorded, but its time ends with its own last one.
    const std::vector<std::int64_t> closedAt{4, 5, 9, 11};
    const std::vector<std::int64_t> sampled{4, 1, 4, 2};
    const std::vector<std::string> wallTimes{"2026-10-11T12:00:02.500Z", "2026-10-11T12:00:03.000Z", "2026-10-11T12:00:05.000Z", "2026-10-11T12:00:05.500Z"};
    const std::vector<std::size_t> rows{1, 3, 4, 6};
    for (std::size_t w = 0; w < 4; ++w) {
        const auto& closed = e[rows[w]];
        EXPECT_EQ(closed.window, static_cast<std::int64_t>(w));
        EXPECT_EQ(closed.iteration, closedAt[w]);
        EXPECT_EQ(closed.phase, w < 2 ? "burn_in" : "sampling");
        EXPECT_EQ(closed.wallTime, wallTimes[w]);
        EXPECT_EQ(closed.payload.at("sampling_iterations"), sampled[w]) << w;
        EXPECT_DOUBLE_EQ(closed.payload.at("seconds"), 0.5 * static_cast<double>(sampled[w])) << w;
        EXPECT_DOUBLE_EQ(closed.payload.at("iterations_per_second"), 2.0) << w;
    }

    EXPECT_EQ(e[5].iteration, 11);
    EXPECT_EQ(e[5].window, 3);
    EXPECT_EQ(e[5].payload, (nlohmann::json{{"reason", "completed"}, {"seconds", 5.5}}));
}

// -t 5, so iterations 3 and 7 write no Trace row. A value that stays non-finite is raised once, when it becomes so.
TEST_F(RunRecorderTest, NonFiniteValueIsRaisedAtTheSamplingIterationItAppearsAtEvenOffTheThinGrid) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    auto recorder    = makeRecorder(5, 10);
    recorder.recordStart(ensemble());
    for (int k = 1; k <= 8; ++k) {
        slots[2].likelihood = (k == 3 or k == 4) ? nan : -30.0;
        slots[0].prior      = k >= 7 ? -inf : -1.0;
        slots[3].likelihood = k == 8 ? inf : -40.0;
        recorder.recordIteration(ensemble());
    }
    recorder.recordEnd();

    std::vector<Event> raised;
    for (const auto& e : events(sink->table("events"))) {
        if (e.type == "non_finite_value") {
            raised.push_back(e);
        }
    }
    ASSERT_EQ(raised.size(), 3u);
    EXPECT_EQ(raised[0].iteration, 3);
    EXPECT_EQ(raised[0].payload, (nlohmann::json{{"temperature_slot", 2}, {"term", "likelihood"}, {"value", "nan"}}));
    EXPECT_EQ(raised[1].iteration, 7);
    EXPECT_EQ(raised[1].payload, (nlohmann::json{{"temperature_slot", 0}, {"term", "prior"}, {"value", "-inf"}}));
    EXPECT_EQ(raised[2].iteration, 8);
    EXPECT_EQ(raised[2].payload, (nlohmann::json{{"temperature_slot", 3}, {"term", "likelihood"}, {"value", "inf"}}));
}

TEST_F(RunRecorderTest, NonFiniteValueIsRaisedForTheStartingState) {
    slots[1].prior = -std::numeric_limits<double>::infinity();
    auto recorder  = makeRecorder(5, 10);
    recorder.recordStart(ensemble());
    recorder.recordEnd();

    const auto e = events(sink->table("events"));
    ASSERT_EQ(types(e), (std::vector<std::string>{"run_started", "non_finite_value", "run_ended", "window_closed"}));
    EXPECT_EQ(e[1].iteration, 0);
    EXPECT_EQ(e[1].payload, (nlohmann::json{{"temperature_slot", 1}, {"term", "prior"}, {"value", "-inf"}}));
}

// -t 1, --window 10. Replica exchange adapts its ladder after burn-in iterations 1 and 2; the second comes out NaN and
// is not applied. Burn-in ends after iteration 2.
TEST_F(RunRecorderTest, TemperatureAdaptationsAndBurnInEndAreRecordedInBurnIn) {
    auto recorder = makeRecorder(1, 10);
    recorder.recordStart(ensemble());
    recorder.recordIteration(ensemble());
    recorder.recordTemperatureAdaptation({.applied = true, .temperatures = {1.0, 0.7, 0.5, 0.2}, .barriers = {0.25, 0.5, 0.75}});
    recorder.recordIteration(ensemble());
    recorder.recordTemperatureAdaptation({.applied = false, .temperatures = {1.0, 0.7, 0.5, 0.2}, .barriers = {0.0, 1.0, 0.5}});
    recorder.recordBurnInEnd(finalAdaptation());
    recorder.recordIteration(ensemble());
    recorder.recordEnd();

    const auto e = events(sink->table("events"));
    ASSERT_EQ(types(e), (std::vector<std::string>{"run_started", "temperatures_adapted", "temperatures_adapted", "burn_in_ended", "window_closed", "run_ended", "window_closed"}));
    EXPECT_EQ(e[1].iteration, 1);
    EXPECT_EQ(e[1].phase, "burn_in");
    EXPECT_EQ(e[1].payload, (nlohmann::json{{"applied", true}, {"temperatures", {1.0, 0.7, 0.5, 0.2}}, {"barriers", {0.25, 0.5, 0.75}}}));
    EXPECT_EQ(e[2].iteration, 2);
    EXPECT_EQ(e[2].payload.at("applied"), false);

    // Proposal variances by sampler id; SALT's are by Allele. A sampler that doesn't adapt is left out.
    EXPECT_EQ(e[3].iteration, 2);
    EXPECT_EQ(e[3].phase, "burn_in");
    EXPECT_EQ(e[3].window, 0);
    EXPECT_EQ(e[3].payload, (nlohmann::json{{"proposal_variances", {{"ConstrainedContinuousRandomWalk/mean_coi", {0.5}}, {"SALTSampler/allele_frequencies/AS1", {0.1, 0.2, 0.3}}}},
                                            {"temperatures", {1.0, 0.7, 0.5, 0.2}}}));
    EXPECT_EQ(e[5].phase, "sampling");
    EXPECT_EQ(e[5].window, 1);
}

// The Run (5 burn-in, -t 2, --window 4) was checkpointed after sampling iteration 7 closed Window 2, and is resumed with
// -s grown to 10. Numbering continues from there: no starting-state row, the `-t` and Window grids of the sampling
// Phase, and Blocks counted from the restored counts.
TEST_F(RunRecorderTest, AResumedRunContinuesIterationAndWindowNumbering) {
    acceptance[0][0] = {.accepted = 10, .proposed = 20};
    auto recorder    = makeRecorder(2, 4);
    recorder.recordResume({.position = {.iteration = 7, .phase = Phase::Sampling, .lastClosedWindow = 2}, .checkpoint = "run/checkpoint.cbor", .samples = 10}, ensemble());
    for (int k = 8; k <= 11; ++k) {
        acceptance[0][0].accepted += 1;
        acceptance[0][0].proposed += 1;
        recorder.recordIteration(ensemble());
    }
    recorder.recordEnd();

    EXPECT_EQ(sink->windows(), (std::vector<WindowMarker>{
                                       {.window = 3, .firstIteration = 8, .lastIteration = 9, .phase = "sampling"},
                                       {.window = 4, .firstIteration = 10, .lastIteration = 11, .phase = "sampling"},
                               }));
    EXPECT_EQ(sink->table("run_scalars").column<std::int64_t>("iteration"), (std::vector<std::int64_t>{9, 11}));
    EXPECT_EQ(sink->table("acceptance").column<std::int64_t>("accepted"), (std::vector<std::int64_t>{2, 0, 0, 2, 0, 0}));

    const auto e = events(sink->table("events"));
    ASSERT_EQ(types(e), (std::vector<std::string>{"run_resumed", "window_closed", "run_ended", "window_closed"}));
    EXPECT_EQ(e[0].iteration, 7);
    EXPECT_EQ(e[0].phase, "sampling");
    EXPECT_EQ(e[0].window, 3);
    EXPECT_EQ(e[0].payload, (nlohmann::json{{"iteration", 7}, {"window", 2}, {"checkpoint", "run/checkpoint.cbor"}, {"samples", 10}}));
}

TEST_F(RunRecorderTest, BlocksHaveNoStartingStateRow) {
    auto recorder = makeRecorder(2, 4);
    recorder.recordStart(ensemble());
    recorder.recordEnd();

    EXPECT_EQ(sink->table("acceptance").rows(), 0u);
    EXPECT_EQ(sink->table("swaps").rows(), 0u);
}

namespace {
    struct Committed {
        RunPosition position;
        bool burnInEnded;

        bool operator==(const Committed&) const = default;
    };

    // Every committed Window covers at least one Sampling iteration.
    void expectNoIterationlessWindow(const std::vector<WindowMarker>& windows) {
        for (const auto& w : windows) {
            EXPECT_LE(w.firstIteration, w.lastIteration) << "Window " << w.window;
        }
    }
}// namespace

// Each committed Window reports where the Run stands, with the Phase of the next Sampling iteration and what remains of
// the Run (5 burn-in and 6 sampling iterations), so a checkpoint can be taken there. The Window burn-in end is
// recorded in says so. The Window a second interrupt closes is not reported: the State may be mid Sampling iteration.
TEST_F(RunRecorderTest, EachCommittedWindowIsReportedExceptOneASecondInterruptCloses) {
    auto recorder = makeRecorder(1, 2);
    std::vector<Committed> committed;
    recorder.onWindowCommitted([&](const WindowCommit& c) { committed.push_back({c.position, c.burnInEnded}); });
    recorder.recordStart(ensemble());
    for (int k = 1; k <= 5; ++k) {
        recorder.recordIteration(ensemble());
    }
    recorder.recordBurnInEnd(finalAdaptation());
    recorder.recordIteration(ensemble());
    recorder.recordInterrupt("SIGINT", Interrupt::First);
    recorder.recordInterrupt("SIGTERM", Interrupt::Second);

    EXPECT_EQ(sink->windows(), (std::vector<WindowMarker>{
                                       {.window = 0, .firstIteration = 0, .lastIteration = 2, .phase = "burn_in"},
                                       {.window = 1, .firstIteration = 3, .lastIteration = 4, .phase = "burn_in"},
                                       {.window = 2, .firstIteration = 5, .lastIteration = 5, .phase = "burn_in"},
                                       {.window = 3, .firstIteration = 6, .lastIteration = 6, .phase = "sampling"},
                               }));
    EXPECT_EQ(committed, (std::vector<Committed>{
                                 {{.iteration = 2, .phase = Phase::BurnIn, .lastClosedWindow = 0, .burnInRemaining = 3, .samplingRemaining = 6}, false},
                                 {{.iteration = 4, .phase = Phase::BurnIn, .lastClosedWindow = 1, .burnInRemaining = 1, .samplingRemaining = 6}, false},
                                 {{.iteration = 5, .phase = Phase::Sampling, .lastClosedWindow = 2, .burnInRemaining = 0, .samplingRemaining = 6}, true},
                         }));
}

// -t 1, --window 2. A Window that fills stays open until the next Sampling iteration is recorded, so a Run event raised
// after its last Sampling iteration goes into it.
TEST_F(RunRecorderTest, AFullWindowStaysOpenUntilTheNextSamplingIteration) {
    auto recorder = makeRecorder(1, 2);
    recorder.recordStart(ensemble());
    recorder.recordIteration(ensemble());
    recorder.recordIteration(ensemble());
    EXPECT_TRUE(recorder.windowComplete());
    recorder.recordTemperatureAdaptation({.applied = true, .temperatures = {1.0, 0.7, 0.5, 0.2}, .barriers = {0.25, 0.5, 0.75}});
    EXPECT_TRUE(sink->windows().empty());

    recorder.recordIteration(ensemble());
    EXPECT_FALSE(recorder.windowComplete());
    ASSERT_EQ(sink->windows(), (std::vector<WindowMarker>{{.window = 0, .firstIteration = 0, .lastIteration = 2, .phase = "burn_in"}}));
    const auto e = events(sink->windowTables(0).at("events"));
    ASSERT_EQ(types(e), (std::vector<std::string>{"run_started", "temperatures_adapted", "window_closed"}));
    EXPECT_EQ(e[1].iteration, 2);
}

// -t 1, --window 2, 4 burn-in and 2 sampling iterations: burn-in and the Run both end on a Window boundary. Their
// events go into the full Window rather than a Window of their own, and the Window burn-in end is recorded in is the
// one a checkpoint at burn-in end covers.
TEST_F(RunRecorderTest, BurnInEndAndRunEndOnAWindowBoundaryGoIntoThatWindow) {
    auto manifest    = exampleManifest();
    manifest.burnIn  = 4;
    manifest.samples = 2;
    auto recorder    = makeRecorder(1, 2, 1e-9, manifest);
    std::vector<Committed> committed;
    recorder.onWindowCommitted([&](const WindowCommit& c) { committed.push_back({c.position, c.burnInEnded}); });
    recorder.recordStart(ensemble());
    for (int k = 1; k <= 4; ++k) {
        recorder.recordIteration(ensemble());
    }
    recorder.recordBurnInEnd(finalAdaptation());
    EXPECT_EQ(sink->windows().size(), 1u);
    for (int k = 5; k <= 6; ++k) {
        recorder.recordIteration(ensemble());
    }
    recorder.recordEnd();

    EXPECT_EQ(sink->windows(), (std::vector<WindowMarker>{
                                       {.window = 0, .firstIteration = 0, .lastIteration = 2, .phase = "burn_in"},
                                       {.window = 1, .firstIteration = 3, .lastIteration = 4, .phase = "burn_in"},
                                       {.window = 2, .firstIteration = 5, .lastIteration = 6, .phase = "sampling"},
                               }));
    expectNoIterationlessWindow(sink->windows());
    EXPECT_EQ(types(events(sink->windowTables(1).at("events"))), (std::vector<std::string>{"burn_in_ended", "window_closed"}));
    EXPECT_EQ(types(events(sink->windowTables(2).at("events"))), (std::vector<std::string>{"run_ended", "window_closed"}));
    ASSERT_EQ(committed.size(), 3u);
    EXPECT_EQ(committed[1], (Committed{{.iteration = 4, .phase = Phase::Sampling, .lastClosedWindow = 1, .burnInRemaining = 0, .samplingRemaining = 2}, true}));
    EXPECT_FALSE(committed[2].burnInEnded);
}

// With no burn-in and no sampling iterations, the starting State's Window holds burn-in end and the Run's end.
TEST_F(RunRecorderTest, ARunWithNoSamplingIterationsCommitsOneWindow) {
    auto recorder = makeRecorder(1, 2);
    recorder.recordStart(ensemble());
    recorder.recordBurnInEnd(finalAdaptation());
    recorder.recordEnd();

    EXPECT_EQ(sink->windows(), (std::vector<WindowMarker>{{.window = 0, .firstIteration = 0, .lastIteration = 0, .phase = "burn_in"}}));
    EXPECT_EQ(types(events(sink->table("events"))), (std::vector<std::string>{"run_started", "burn_in_ended", "run_ended", "window_closed"}));
}

// A resumed Run that ends before recording a Sampling iteration commits nothing: the store stays as its checkpoint
// left it.
TEST_F(RunRecorderTest, AResumedRunThatRecordsNoSamplingIterationCommitsNothing) {
    auto recorder = makeRecorder(2, 4);
    std::vector<Committed> committed;
    recorder.onWindowCommitted([&](const WindowCommit& c) { committed.push_back({c.position, c.burnInEnded}); });
    recorder.recordResume({.position = {.iteration = 7, .phase = Phase::Sampling, .lastClosedWindow = 2}, .checkpoint = "run/checkpoint.cbor", .samples = 6}, ensemble());
    recorder.recordInterrupt("SIGINT", Interrupt::First);
    recorder.recordEnd();

    EXPECT_TRUE(sink->windows().empty());
    EXPECT_FALSE(sink->hasUncommittedTables());
    EXPECT_TRUE(committed.empty());
}

// The Parent-set ids a resumed Run is handed are those its checkpoint carried: a Parent set seen before the checkpoint
// keeps its id and is not listed again.
TEST_F(RunRecorderTest, AResumedRunKeepsTheParentSetIdsItIsHanded) {
    replica.makeEveryParentSetPossible();
    ParentSetIdDictionary carried;
    {
        auto before = makeRecorder(1, 10, 0.0);
        before.onWindowCommitted([&](const WindowCommit& c) { carried = c.parentSetIds; });
        before.recordStart(ensemble());
        before.recordEnd();
    }
    ASSERT_EQ(carried.at("3").ids.size(), 7u);
    EXPECT_EQ(carried.at("3").nextId, 7u);

    auto recorder = makeRecorder(1, 10, 0.0);
    recorder.recordResume({.position = {.iteration = 7, .phase = Phase::Sampling, .lastClosedWindow = 2}, .checkpoint = "run/checkpoint.cbor", .samples = 10, .parentSetIds = carried}, ensemble());
    recorder.recordIteration(ensemble());
    recorder.recordEnd();

    EXPECT_EQ(sink->table("parent_set_members").rows(), 0u);
    const auto sets = sink->table("parent_sets");
    ASSERT_GT(sets.rows(), 0u);
    std::set<std::pair<std::string, std::int64_t>> known;
    for (const auto& [infection, ids] : carried) {
        for (const auto& [members, id] : ids.ids) {
            known.insert({infection, static_cast<std::int64_t>(id)});
        }
    }
    for (std::size_t r = 0; r < sets.rows(); ++r) {
        EXPECT_TRUE(known.contains({sets.column<std::string>("infection")[r], sets.column<std::int64_t>("set_id")[r]}));
    }
}
