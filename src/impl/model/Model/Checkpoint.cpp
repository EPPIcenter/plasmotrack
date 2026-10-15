#include "Checkpoint.h"

#include "core/io/DurableFile.h"
#include "core/io/OutputFailure.h"

#include <fmt/core.h>

#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <bit>

#include <cerrno>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <system_error>

namespace transmission_nets::impl::Model::checkpoint {
    namespace fs = std::filesystem;
    using nlohmann::json;

    namespace {
        // A checkpoint that can't be used is a refusal, which exits 1.
        [[noreturn]] void refuse(std::string message) {
            throw Refused({.kind = RefusalKind::UnusableCheckpoint, .message = std::move(message)});
        }

        // The saved value for `key`; a missing key means the checkpoint was made from another input.
        template<typename Value>
        const Value& savedFor(const std::map<std::string, Value>& saved, const std::string& key, const std::string& what) {
            const auto it = saved.find(key);
            if (it == saved.end()) {
                refuse(fmt::format("checkpoint has no {} {}; it was made from another input", what, key));
            }
            return it->second;
        }

        void requireSameNames(const std::vector<std::string>& saved, const std::vector<std::string>& current, const std::string& what) {
            if (saved != current) {
                refuse(fmt::format("checkpoint has other {} than this input; it was made from another input", what));
            }
        }

        template<typename Value>
        void requireLength(const std::vector<Value>& saved, const size_t expected, const std::string& what) {
            if (saved.size() != expected) {
                refuse(fmt::format("checkpoint has {} {}, not {}; it was made from another input", saved.size(), what, expected));
            }
        }

        std::vector<std::string> infectionIds(const std::vector<std::shared_ptr<InfectionEvent>>& infections) {
            std::vector<std::string> ids;
            ids.reserve(infections.size());
            for (const auto& infection : infections) {
                ids.push_back(infection->id());
            }
            return ids;
        }

        std::vector<std::string> locusLabels(const State& state) {
            std::vector<std::string> labels;
            for (const auto& [label, locus] : state.loci) {
                labels.push_back(label);
            }
            return labels;
        }

        // Appends Allele bits to bytes, filling each byte from its most significant bit.
        class BitWriter {
        public:
            void push(const bool bit) {
                if (bits_ % 8 == 0) {
                    bytes_.push_back(0);
                }
                bytes_.back() |= static_cast<std::uint8_t>(static_cast<unsigned>(bit) << (7 - bits_ % 8));
                ++bits_;
            }
            std::vector<std::uint8_t> take() { return std::move(bytes_); }

        private:
            std::vector<std::uint8_t> bytes_;
            size_t bits_ = 0;
        };

        class BitReader {
        public:
            explicit BitReader(const std::vector<std::uint8_t>& bytes) : bytes_(bytes) {}
            bool next() {
                if (bits_ / 8 >= bytes_.size()) {
                    refuse("checkpoint has fewer Allele bits than this input's genotypes; it was made from another input");
                }
                const bool bit = (bytes_[bits_ / 8] >> (7 - bits_ % 8)) & 1;
                ++bits_;
                return bit;
            }
            void requireAllRead() const {
                if ((bits_ + 7) / 8 != bytes_.size()) {
                    refuse("checkpoint has more Allele bits than this input's genotypes; it was made from another input");
                }
            }

        private:
            const std::vector<std::uint8_t>& bytes_;
            size_t bits_ = 0;
        };

        // Visits a State's Quantities into a ChainState.
        struct Capturing {
            const State& state;
            ChainState saved{};
            BitWriter genotypes{};

            void scalar(const char* name, const std::shared_ptr<core::parameters::Parameter<double>>& p) {
                saved.scalars.emplace(name, p->value());
            }
            void perInfection(const char* name, size_t, const std::shared_ptr<core::parameters::Parameter<double>>& p) {
                saved.perInfection[name].push_back(p->value());
            }
            void alleleFrequencies(const std::shared_ptr<LocusImpl>&, const std::shared_ptr<core::parameters::Parameter<AlleleFrequencyImpl>>& p) {
                saved.alleleFrequencies.push_back(p->value().frequencies());
            }
            void latentGenotype(size_t, const std::shared_ptr<LocusImpl>&, const std::shared_ptr<core::parameters::Parameter<GeneticsImpl>>& p) {
                const auto& alleles = p->value();
                for (size_t allele = 0; allele < alleles.totalAlleles(); ++allele) {
                    genotypes.push(alleles.allele(allele));
                }
            }
            void ordering(const std::shared_ptr<OrderingImpl>& ordering) {
                std::map<const InfectionEvent*, std::uint32_t> position;
                for (std::uint32_t ii = 0; ii < state.infections.size(); ++ii) {
                    position.emplace(state.infections[ii].get(), ii);
                }
                for (const auto& infection : ordering->value()) {
                    saved.ordering.push_back(position.at(infection.get()));
                }
            }
            ChainState take() {
                saved.latentGenotypes = genotypes.take();
                return std::move(saved);
            }
        };

        // Visits a ChainState into a State no Model has been built on yet.
        struct Restoring {
            State& state;
            const ChainState& saved;
            BitReader genotypes;
            size_t locusAt = 0;

            void scalar(const char* name, const std::shared_ptr<core::parameters::Parameter<double>>& p) {
                p->initializeValue(savedFor(saved.scalars, name, "Quantity"));
            }
            void perInfection(const char* name, const size_t index, const std::shared_ptr<core::parameters::Parameter<double>>& p) {
                const auto& values = savedFor(saved.perInfection, name, "per-Infection Quantity");
                if (index == 0) {
                    requireLength(values, state.infections.size(), name);
                }
                p->initializeValue(values[index]);
            }
            void alleleFrequencies(const std::shared_ptr<LocusImpl>&, const std::shared_ptr<core::parameters::Parameter<AlleleFrequencyImpl>>& p) {
                p->initializeValue(AlleleFrequencyImpl::exactly(saved.alleleFrequencies[locusAt++]));
            }
            void latentGenotype(size_t, const std::shared_ptr<LocusImpl>&, const std::shared_ptr<core::parameters::Parameter<GeneticsImpl>>& p) {
                GeneticsImpl alleles(static_cast<int>(p->value().totalAlleles()));
                for (size_t allele = 0; allele < alleles.totalAlleles(); ++allele) {
                    alleles.set(allele, genotypes.next());
                }
                p->initializeValue(alleles);
            }
            void ordering(const std::shared_ptr<OrderingImpl>& ordering) {
                genotypes.requireAllRead();
                std::vector<std::shared_ptr<InfectionEvent>> order;
                for (const auto at : saved.ordering) {
                    if (at >= state.infections.size()) {
                        refuse(fmt::format("checkpoint's Ordering has Infection {} of {}", at, state.infections.size()));
                    }
                    order.push_back(state.infections[at]);
                }
                ordering->restoreOrder(std::move(order));
            }
        };

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

        json toJson(const Ensemble::ExchangeState& exchange) {
            json swapCounts = json::array();
            for (const auto& count : exchange.swapCounts) {
                swapCounts.push_back({{"accepted", count.accepted}, {"attempted", count.attempted}});
            }
            return {{"replicas", exchange.replicas},
                    {"temperatures", exchange.temperatures},
                    {"swap_barriers", exchange.swapBarriers},
                    {"barrier_iterations", exchange.barrierIterations},
                    {"second_pair_set_next", exchange.secondPairSetNext},
                    {"swap_counts", swapCounts}};
        }

        Ensemble::ExchangeState exchangeFromJson(const json& j) {
            Ensemble::ExchangeState exchange{.replicas          = j.at("replicas"),
                                             .temperatures      = j.at("temperatures"),
                                             .swapBarriers      = j.at("swap_barriers"),
                                             .barrierIterations = j.at("barrier_iterations"),
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
                refuse("checkpoint has sampler counts of unequal lengths");
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
                refuse(fmt::format("checkpoint has a replica at scheduler iteration {}, out of range", iterations));
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
            refuse("checkpoint has an unknown Phase " + name);
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

        // Writes CBOR straight to a file through a small buffer, so a checkpoint is never held whole in memory. The
        // small parts go through nlohmann_json; the large arrays are written item by item.
        class CborFile {
        public:
            explicit CborFile(fs::path path) : path_(std::move(path)), fd_(::open(path_.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644)) {
                if (fd_ < 0) {
                    core::io::failOutput("open", path_, std::strerror(errno));
                }
                buffer_.reserve(BUFFER_BYTES);
            }
            CborFile(const CborFile&)            = delete;
            CborFile& operator=(const CborFile&) = delete;
            ~CborFile() {
                if (fd_ >= 0) {
                    ::close(fd_);
                }
            }

            void map(const size_t entries) { head(5, entries); }
            void array(const size_t items) { head(4, items); }
            void text(const std::string_view value) {
                head(3, value.size());
                append(value);
            }
            void bytes(const std::vector<std::uint8_t>& value) {
                head(2, value.size());
                append({reinterpret_cast<const char*>(value.data()), value.size()});
            }
            void number(const std::uint64_t value) { head(0, value); }
            void number(const std::int64_t value) {
                if (value < 0) {
                    head(1, static_cast<std::uint64_t>(-(value + 1)));
                } else {
                    head(0, static_cast<std::uint64_t>(value));
                }
            }
            void number(const double value) {
                const auto bits = std::bit_cast<std::uint64_t>(value);
                buffer_.push_back(static_cast<char>(0xfb));
                for (int shift = 56; shift >= 0; shift -= 8) {
                    buffer_.push_back(static_cast<char>((bits >> shift) & 0xff));
                }
                flushIfFull();
            }
            template<typename Number>
            void numbers(const std::vector<Number>& values) {
                array(values.size());
                for (const auto value : values) {
                    number(value);
                }
            }
            void texts(const std::vector<std::string>& values) {
                array(values.size());
                for (const auto& value : values) {
                    text(value);
                }
            }
            void value(const json& small) {
                const auto encoded = json::to_cbor(small);
                append({reinterpret_cast<const char*>(encoded.data()), encoded.size()});
            }

            void close() {
                flush();
                const int fd = fd_;
                fd_          = -1;
                if (::close(fd) != 0) {
                    core::io::failOutput("write", path_, std::strerror(errno));
                }
            }

        private:
            static constexpr size_t BUFFER_BYTES = size_t{1} << 20;

            void head(const std::uint8_t major, const std::uint64_t value) {
                const auto type = static_cast<std::uint8_t>(major << 5);
                if (value < 24) {
                    buffer_.push_back(static_cast<char>(type | value));
                } else {
                    const int width = value <= 0xff ? 1 : value <= 0xffff ? 2 : value <= 0xffffffff ? 4 : 8;
                    buffer_.push_back(static_cast<char>(type | (width == 1 ? 24 : width == 2 ? 25 : width == 4 ? 26 : 27)));
                    for (int shift = 8 * (width - 1); shift >= 0; shift -= 8) {
                        buffer_.push_back(static_cast<char>((value >> shift) & 0xff));
                    }
                }
                flushIfFull();
            }

            void append(const std::string_view raw) {
                buffer_.append(raw);
                flushIfFull();
            }

            void flushIfFull() {
                if (buffer_.size() >= BUFFER_BYTES) {
                    flush();
                }
            }

            void flush() {
                size_t written = 0;
                while (written < buffer_.size()) {
                    const ssize_t n = ::write(fd_, buffer_.data() + written, buffer_.size() - written);
                    if (n < 0) {
                        if (errno == EINTR) {
                            continue;
                        }
                        core::io::failOutput("write", path_, std::strerror(errno));
                    }
                    written += static_cast<size_t>(n);
                }
                buffer_.clear();
            }

            fs::path path_;
            int fd_;
            std::string buffer_;
        };

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

    Layout stateLayout(const State& state) {
        return {.infections = infectionIds(state.infections), .latentParents = infectionIds(state.latentParents), .loci = locusLabels(state)};
    }

    ChainState captureState(const State& state) {
        Capturing capturing{.state = state};
        state.quantities(capturing);
        return capturing.take();
    }

    // Fills a State no Model has been built on yet, from the same input as `layout`.
    void restoreState(const Layout& layout, const ChainState& saved, State& state) {
        requireSameNames(layout.infections, infectionIds(state.infections), "Infections");
        requireSameNames(layout.latentParents, infectionIds(state.latentParents), "Latent parents");
        requireSameNames(layout.loci, locusLabels(state), "Loci");
        requireLength(saved.alleleFrequencies, state.loci.size(), "Allele frequencies");
        requireLength(saved.ordering, state.infections.size(), "Infections in the Ordering");
        Restoring restoring{.state = state, .saved = saved, .genotypes = BitReader(saved.latentGenotypes)};
        state.quantities(restoring);
    }


    Checkpoint capture(const Ensemble& ensemble, Identity identity, RunPosition position, ParentSetIdDictionary parentSetIds) {
        Checkpoint checkpoint{.identity = std::move(identity), .position = position, .parentSetIds = std::move(parentSetIds), .layout = stateLayout(*ensemble.chains.front().state)};
        checkpoint.ensemble = ensemble.record([&](const State& state) {
            if (infectionIds(state.infections) != checkpoint.layout.infections) {
                throw std::logic_error("chains of one ensemble have different Infections");
            }
            return captureState(state);
        });
        return checkpoint;
    }

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
            refuse(fmt::format("cannot open the checkpoint {}: {}", file.string(), std::strerror(errno)));
        }
        const std::vector<std::uint8_t> bytes{std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
        try {
            const json saved = json::from_cbor(bytes);
            const int version = saved.at("format_version");
            if (version != FORMAT_VERSION) {
                refuse(fmt::format("{} has checkpoint format version {}; this engine reads version {}", file.string(), version, FORMAT_VERSION));
            }
            return fromJson(saved);
        } catch (const json::exception& e) {
            refuse(fmt::format("{} is not a readable checkpoint: {}", file.string(), e.what()));
        }
    }

    std::unique_ptr<Ensemble> restore(const Checkpoint& checkpoint,
                                      const std::shared_ptr<Engine>& runRng,
                                      const nlohmann::json& input,
                                      const std::vector<core::computation::Probability>& symptomaticIDPDist,
                                      const std::vector<core::computation::Probability>& asymptomaticIDPDist,
                                      const bool nullModel,
                                      const unsigned int numCores) {
        std::unique_ptr<Ensemble> ensemble;
        try {
            ensemble = std::make_unique<Ensemble>(checkpoint.ensemble, runRng, numCores, [&](const size_t, const ChainState& saved, std::shared_ptr<Engine> chainRng) {
                auto state = std::make_shared<State>(input, symptomaticIDPDist, asymptomaticIDPDist, std::move(chainRng), nullModel);
                restoreState(checkpoint.layout, saved, *state);
                return state;
            });
        } catch (const std::invalid_argument& e) {
            refuse(fmt::format("checkpoint of Run {} does not fit this input: {}", checkpoint.identity.runId, e.what()));
        }

        // The burn-in iteration count isn't saved: in burn-in it is the position's iteration. Once burn-in has ended
        // it is never read.
        auto exchange             = ensemble->exchangeState();
        exchange.burnInIterations = static_cast<int>(checkpoint.position.iteration);
        // Once burn-in has ended, barriers still pending can't be told apart by pair: the ladder works out each pair's
        // attempts from burn-in iterations that run back to back, and Sampling iterations have alternated the pairs
        // since. A Run started from here begins its own.
        if (checkpoint.position.phase == Phase::Sampling) {
            std::ranges::fill(exchange.swapBarriers, 0.0);
            exchange.barrierIterations = 0;
        }
        ensemble->restoreExchangeState(exchange);
        return ensemble;
    }
}// namespace transmission_nets::impl::Model::checkpoint
