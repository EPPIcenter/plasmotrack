//
// The Checkpoint codec: a Checkpoint written as CBOR and read back. See Checkpoint.h.
//

#include "Checkpoint.h"

#include "core/io/CborFile.h"
#include "core/io/DurableFile.h"

#include <fmt/core.h>

#include <cerrno>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <stdexcept>

namespace transmission_nets::impl::Model::checkpoint {
    namespace fs = std::filesystem;
    using nlohmann::json;
    using core::io::CborFile;

    void refuseCheckpoint(std::string message) {
        throw Refused({.kind = RefusalKind::UnusableCheckpoint, .message = std::move(message)});
    }

    namespace {
        std::vector<std::uint8_t> bytesFrom(const json& j) {
            const auto& bytes = j.get_binary();
            return {bytes.begin(), bytes.end()};
        }

        ChainState chainStateFromJson(const json& j) {
            return {.scalars           = j.at("scalars"),
                    .perInfection      = j.at("per_infection"),
                    .alleleFrequencies = j.at("allele_frequencies"),
                    .latentGenotypes   = bytesFrom(j.at("latent_genotypes")),
                    .ordering          = j.at("ordering")};
        }

        json toJson(const core::samplers::ExchangeState& exchange) {
            json swapCounts = json::array();
            for (const auto& count : exchange.swapCounts) {
                swapCounts.push_back({{"accepted", count.accepted}, {"attempted", count.attempted}});
            }
            return {{"replicas", exchange.replicas},
                    {"temperatures", exchange.temperatures},
                    {"swap_barriers", exchange.swapBarriers},
                    {"barrier_iterations", exchange.barrierIterations},
                    {"burn_in_iterations", exchange.burnInIterations},
                    {"second_pair_set_next", exchange.secondPairSetNext},
                    {"swap_counts", swapCounts}};
        }

        core::samplers::ExchangeState exchangeFromJson(const json& j) {
            core::samplers::ExchangeState exchange{.replicas          = j.at("replicas"),
                                             .temperatures      = j.at("temperatures"),
                                             .swapBarriers      = j.at("swap_barriers"),
                                             .barrierIterations = j.at("barrier_iterations"),
                                             .burnInIterations  = j.at("burn_in_iterations"),
                                             .secondPairSetNext = j.at("second_pair_set_next")};
            for (const auto& count : j.at("swap_counts")) {
                exchange.swapCounts.push_back({.accepted = count.at("accepted"), .attempted = count.at("attempted")});
            }
            return exchange;
        }

        std::vector<core::samplers::SamplerCounts> countsFromJson(const json& j) {
            const auto accepted = j.at("accepted").get<std::vector<std::uint64_t>>();
            const auto proposed = j.at("proposed").get<std::vector<std::uint64_t>>();
            if (accepted.size() != proposed.size()) {
                refuseCheckpoint("checkpoint has sampler counts of unequal lengths");
            }
            std::vector<core::samplers::SamplerCounts> counts;
            counts.reserve(accepted.size());
            for (size_t ii = 0; ii < accepted.size(); ++ii) {
                counts.push_back({.accepted = accepted[ii], .proposed = proposed[ii]});
            }
            return counts;
        }

        std::map<std::uint32_t, core::samplers::AdaptationState> adaptationFromJson(const json& j) {
            const auto& samplers  = j.at("samplers");
            const auto& updates   = j.at("updates");
            const auto& variances = j.at("variances");
            const auto& accepted  = j.at("accepted");
            const auto& proposed  = j.at("proposed");
            std::map<std::uint32_t, core::samplers::AdaptationState> adaptation;
            for (size_t ii = 0; ii < samplers.size(); ++ii) {
                auto& state   = adaptation[samplers.at(ii).get<std::uint32_t>()];
                state.updates = updates.at(ii);
                for (size_t kk = 0; kk < variances.at(ii).size(); ++kk) {
                    state.scales.push_back({.variance = variances.at(ii).at(kk), .counts = {.accepted = accepted.at(ii).at(kk), .proposed = proposed.at(ii).at(kk)}});
                }
            }
            return adaptation;
        }

        Replica replicaFromJson(const json& j) {
            // The Sampler set counts in int.
            const auto iterations = j.at("scheduler_iterations").get<std::int64_t>();
            if (iterations < 0 or iterations > std::numeric_limits<int>::max()) {
                refuseCheckpoint(fmt::format("checkpoint has a replica at scheduler iteration {}, out of range", iterations));
            }
            return {.state    = chainStateFromJson(j.at("state")),
                    .rng      = j.at("rng"),
                    .samplers = {.iterations = static_cast<int>(iterations), .counts = countsFromJson(j.at("sampler_counts")), .adaptation = adaptationFromJson(j.at("adaptation"))}};
        }

        json toJson(const Identity& identity) {
            return {{"run_id", identity.runId}, {"input_hash", identity.inputHash}, {"engine_version", identity.engineVersion}};
        }

        Identity identityFromJson(const json& j) {
            return {.runId = j.at("run_id"), .inputHash = j.at("input_hash"), .engineVersion = j.at("engine_version")};
        }

        std::string phaseName(const Phase phase) {
            return phase == Phase::BurnIn ? "burn_in" : "sampling";
        }

        Phase phaseFromName(const std::string& name) {
            if (name == "burn_in") {
                return Phase::BurnIn;
            }
            if (name == "sampling") {
                return Phase::Sampling;
            }
            refuseCheckpoint("checkpoint has an unknown Phase " + name);
        }

        json toJson(const RunPosition& position) {
            return {{"iteration", position.iteration},
                    {"phase", phaseName(position.phase)},
                    {"last_closed_window", position.lastClosedWindow},
                    {"burn_in_remaining", position.burnInRemaining},
                    {"sampling_remaining", position.samplingRemaining}};
        }

        // Every field is required: a checkpoint is taken only once a Window is committed.
        RunPosition positionFromJson(const json& j) {
            return {.iteration         = j.at("iteration").get<std::int64_t>(),
                    .phase             = phaseFromName(j.at("phase")),
                    .lastClosedWindow  = j.at("last_closed_window").get<std::int64_t>(),
                    .burnInRemaining   = j.at("burn_in_remaining").get<std::int64_t>(),
                    .samplingRemaining = j.at("sampling_remaining").get<std::int64_t>()};
        }

        json toJson(const ParentSetIdDictionary& dictionary) {
            json byInfection = json::object();
            for (const auto& [infection, ids] : dictionary) {
                json sets = json::array();
                for (const auto& [members, setId] : ids.ids) {
                    sets.push_back({{"members", members}, {"set_id", setId}});
                }
                byInfection[infection] = {{"sets", sets}, {"next_id", ids.nextId}};
            }
            return byInfection;
        }

        ParentSetIdDictionary parentSetIdsFromJson(const json& j) {
            ParentSetIdDictionary dictionary;
            for (const auto& [infection, ids] : j.items()) {
                auto& entry = dictionary[infection];
                for (const auto& set : ids.at("sets")) {
                    entry.ids[set.at("members").get<std::vector<std::string>>()] = set.at("set_id");
                }
                entry.nextId = ids.at("next_id");
            }
            return dictionary;
        }

        void write(CborFile& out, const Layout& layout, const std::vector<std::string>& samplerIds) {
            out.map(4);
            out.text("infections");
            out.texts(layout.infections);
            out.text("latent_parents");
            out.texts(layout.latentParents);
            out.text("loci");
            out.texts(layout.loci);
            out.text("samplers");
            out.texts(samplerIds);
        }

        void write(CborFile& out, const ChainState& state) {
            out.map(5);
            out.text("scalars");
            out.map(state.scalars.size());
            for (const auto& [name, value] : state.scalars) {
                out.text(name);
                out.number(value);
            }
            out.text("per_infection");
            out.map(state.perInfection.size());
            for (const auto& [name, values] : state.perInfection) {
                out.text(name);
                out.numbers(values);
            }
            out.text("allele_frequencies");
            out.array(state.alleleFrequencies.size());
            for (const auto& frequencies : state.alleleFrequencies) {
                out.numbers(frequencies);
            }
            out.text("latent_genotypes");
            out.bytes(state.latentGenotypes);
            out.text("ordering");
            out.array(state.ordering.size());
            for (const auto at : state.ordering) {
                out.number(std::uint64_t{at});
            }
        }

        // Each sampler's counts, as two arrays in Layout order.
        void write(CborFile& out, const std::vector<core::samplers::SamplerCounts>& counts) {
            out.map(2);
            out.text("accepted");
            out.array(counts.size());
            for (const auto& count : counts) {
                out.number(count.accepted);
            }
            out.text("proposed");
            out.array(counts.size());
            for (const auto& count : counts) {
                out.number(count.proposed);
            }
        }

        // Parallel arrays over the samplers that adapt: which sampler, its update count, and per proposal scale its
        // variance and counts.
        void write(CborFile& out, const std::map<std::uint32_t, core::samplers::AdaptationState>& adaptation) {
            out.map(5);
            out.text("samplers");
            out.array(adaptation.size());
            for (const auto& [sampler, state] : adaptation) {
                out.number(std::uint64_t{sampler});
            }
            out.text("updates");
            out.array(adaptation.size());
            for (const auto& [sampler, state] : adaptation) {
                out.number(state.updates);
            }
            const auto scales = [&](const char* name, const auto& field) {
                out.text(name);
                out.array(adaptation.size());
                for (const auto& [sampler, state] : adaptation) {
                    out.array(state.scales.size());
                    for (const auto& scale : state.scales) {
                        out.number(field(scale));
                    }
                }
            };
            scales("variances", [](const core::samplers::ProposalScale& scale) { return scale.variance; });
            scales("accepted", [](const core::samplers::ProposalScale& scale) { return scale.counts.accepted; });
            scales("proposed", [](const core::samplers::ProposalScale& scale) { return scale.counts.proposed; });
        }

        void write(CborFile& out, const Replica& replica) {
            out.map(5);
            out.text("state");
            write(out, replica.state);
            out.text("rng");
            out.text(replica.rng);
            out.text("sampler_counts");
            write(out, replica.samplers.counts);
            out.text("adaptation");
            write(out, replica.samplers.adaptation);
            out.text("scheduler_iterations");
            out.number(std::int64_t{replica.samplers.iterations});
        }

        void write(CborFile& out, const Checkpoint& checkpoint) {
            out.map(8);
            out.text("format_version");
            out.number(std::int64_t{FORMAT_VERSION});
            out.text("identity");
            out.value(toJson(checkpoint.identity));
            out.text("run_position");
            out.value(toJson(checkpoint.position));
            out.text("parent_set_ids");
            out.value(toJson(checkpoint.parentSetIds));
            out.text("layout");
            write(out, checkpoint.layout, checkpoint.ensemble.samplerIds);
            out.text("run_rng");
            out.text(checkpoint.ensemble.runRng);
            out.text("replicas");
            out.array(checkpoint.ensemble.replicas.size());
            for (const auto& replica : checkpoint.ensemble.replicas) {
                write(out, replica);
            }
            out.text("exchange");
            out.value(toJson(checkpoint.ensemble.exchange));
        }

        Checkpoint fromJson(const json& j) {
            Checkpoint checkpoint{.identity     = identityFromJson(j.at("identity")),
                                  .position     = positionFromJson(j.at("run_position")),
                                  .parentSetIds = parentSetIdsFromJson(j.at("parent_set_ids")),
                                  .layout       = {.infections    = j.at("layout").at("infections"),
                                                   .latentParents = j.at("layout").at("latent_parents"),
                                                   .loci          = j.at("layout").at("loci")}};
            checkpoint.ensemble.samplerIds = j.at("layout").at("samplers").get<std::vector<std::string>>();
            checkpoint.ensemble.runRng     = j.at("run_rng").get<std::string>();
            for (const auto& replica : j.at("replicas")) {
                checkpoint.ensemble.replicas.push_back(replicaFromJson(replica));
            }
            checkpoint.ensemble.exchange = exchangeFromJson(j.at("exchange"));
            return checkpoint;
        }

    }// namespace

    void save(const Checkpoint& checkpoint, const fs::path& file) {
        core::io::writeDurably(file, [&](const fs::path& temporary) {
            CborFile out(temporary);
            write(out, checkpoint);
            out.close();
        });
    }

    Checkpoint load(const fs::path& file) {
        std::ifstream in(file, std::ios::binary);
        if (!in) {
            refuseCheckpoint(fmt::format("cannot open the checkpoint {}: {}", file.string(), std::strerror(errno)));
        }
        const std::vector<std::uint8_t> bytes{std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
        try {
            const json saved = json::from_cbor(bytes);
            const int version = saved.at("format_version");
            if (version != FORMAT_VERSION) {
                refuseCheckpoint(fmt::format("{} has checkpoint format version {}; this engine reads version {}", file.string(), version, FORMAT_VERSION));
            }
            return fromJson(saved);
        } catch (const json::exception& e) {
            refuseCheckpoint(fmt::format("{} is not a readable checkpoint: {}", file.string(), e.what()));
        }
    }

}// namespace transmission_nets::impl::Model::checkpoint
