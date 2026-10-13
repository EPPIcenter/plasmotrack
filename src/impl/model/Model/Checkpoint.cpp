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

        std::string engineState(const Engine& engine) {
            std::ostringstream out;
            out << engine;
            return out.str();
        }

        void setEngineState(Engine& engine, const std::string& state) {
            std::istringstream in(state);
            in >> engine;
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

        void packGenotypes(const InfectionEvent& owner, const State& state, BitWriter& bits) {
            const auto& genotypes = owner.latentGenotype();
            for (const auto& [label, locus] : state.loci) {
                if (const auto it = genotypes.find(locus); it != genotypes.end()) {
                    const auto& alleles = it->second->value();
                    for (size_t allele = 0; allele < alleles.totalAlleles(); ++allele) {
                        bits.push(alleles.allele(allele));
                    }
                }
            }
        }

        void unpackGenotypes(InfectionEvent& owner, const State& state, BitReader& bits) {
            auto& genotypes = owner.latentGenotype();
            for (const auto& [label, locus] : state.loci) {
                if (const auto it = genotypes.find(locus); it != genotypes.end()) {
                    GeneticsImpl alleles(static_cast<int>(it->second->value().totalAlleles()));
                    for (size_t allele = 0; allele < alleles.totalAlleles(); ++allele) {
                        alleles.set(allele, bits.next());
                    }
                    it->second->initializeValue(alleles);
                }
            }
        }

        Layout stateLayout(const State& state) {
            return {.infections = infectionIds(state.infections), .latentParents = infectionIds(state.latentParents), .loci = locusLabels(state)};
        }

        ChainState captureState(const State& state) {
            ChainState saved;
            BitWriter genotypes;
            for (size_t ii = 0; ii < state.infections.size(); ++ii) {
                const auto& infection = *state.infections[ii];
                packGenotypes(infection, state, genotypes);
                saved.infectionDurations.push_back(infection.infectionDuration()->value());
                saved.expectedFalsePositives.push_back(state.expectedFalsePositives[ii]->value());
                saved.expectedFalseNegatives.push_back(state.expectedFalseNegatives[ii]->value());
            }
            saved.latentGenotypes = genotypes.take();
            for (const auto& [label, locus] : state.loci) {
                saved.alleleFrequencies.push_back(state.alleleFrequencies->alleleFrequencies(locus)->value().frequencies());
            }
            std::map<const InfectionEvent*, std::uint32_t> position;
            for (std::uint32_t ii = 0; ii < state.infections.size(); ++ii) {
                position.emplace(state.infections[ii].get(), ii);
            }
            for (const auto& infection : state.infectionEventOrdering->value()) {
                saved.ordering.push_back(position.at(infection.get()));
            }
            saved.meanCoi                = state.meanCOI->value();
            saved.meanStrainsTransmitted = state.meanStrainsTransmitted->value();
            saved.parentSetSizeProb      = state.parentSetSizeProb->value();
            return saved;
        }

        // Fills a State no Model has been built on yet, from the same input as `layout`.
        void restoreState(const Layout& layout, const ChainState& saved, State& state) {
            requireSameNames(layout.infections, infectionIds(state.infections), "Infections");
            requireSameNames(layout.latentParents, infectionIds(state.latentParents), "Latent parents");
            requireSameNames(layout.loci, locusLabels(state), "Loci");
            const auto numInfections = state.infections.size();
            requireLength(saved.infectionDurations, numInfections, "infection durations");
            requireLength(saved.expectedFalsePositives, numInfections, "expected false positives");
            requireLength(saved.expectedFalseNegatives, numInfections, "expected false negatives");
            requireLength(saved.alleleFrequencies, state.loci.size(), "Allele frequencies");
            requireLength(saved.ordering, numInfections, "Infections in the Ordering");

            BitReader genotypes(saved.latentGenotypes);
            for (size_t ii = 0; ii < numInfections; ++ii) {
                auto& infection = *state.infections[ii];
                unpackGenotypes(infection, state, genotypes);
                infection.infectionDuration()->initializeValue(saved.infectionDurations[ii]);
                state.expectedFalsePositives[ii]->initializeValue(saved.expectedFalsePositives[ii]);
                state.expectedFalseNegatives[ii]->initializeValue(saved.expectedFalseNegatives[ii]);
            }
            genotypes.requireAllRead();
            size_t locusAt = 0;
            for (const auto& [label, locus] : state.loci) {
                state.alleleFrequencies->alleleFrequencies(locus)->initializeValue(AlleleFrequencyImpl::exactly(saved.alleleFrequencies[locusAt++]));
            }
            std::vector<std::shared_ptr<InfectionEvent>> order;
            for (const auto at : saved.ordering) {
                if (at >= numInfections) {
                    refuse(fmt::format("checkpoint's Ordering has Infection {} of {}", at, numInfections));
                }
                order.push_back(state.infections[at]);
            }
            state.infectionEventOrdering->restoreOrder(std::move(order));
            state.meanCOI->initializeValue(saved.meanCoi);
            state.meanStrainsTransmitted->initializeValue(saved.meanStrainsTransmitted);
            state.parentSetSizeProb->initializeValue(saved.parentSetSizeProb);
        }

        std::vector<std::uint8_t> bytesFrom(const json& j) {
            const auto& bytes = j.get_binary();
            return {bytes.begin(), bytes.end()};
        }

        ChainState chainStateFromJson(const json& j) {
            return {.latentGenotypes        = bytesFrom(j.at("latent_genotypes")),
                    .alleleFrequencies      = j.at("allele_frequencies"),
                    .ordering               = j.at("ordering"),
                    .infectionDurations     = j.at("infection_durations"),
                    .expectedFalsePositives = j.at("expected_false_positives"),
                    .expectedFalseNegatives = j.at("expected_false_negatives"),
                    .meanCoi                = j.at("mean_coi"),
                    .meanStrainsTransmitted = j.at("mean_strains_tx"),
                    .parentSetSizeProb      = j.at("parent_set_size_prob")};
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
            return {.state               = chainStateFromJson(j.at("state")),
                    .rng                 = j.at("rng"),
                    .samplerCounts       = countsFromJson(j.at("sampler_counts")),
                    .adaptation          = adaptationFromJson(j.at("adaptation")),
                    .schedulerIterations = j.at("scheduler_iterations")};
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

        void write(CborFile& out, const Layout& layout) {
            out.map(4);
            out.text("infections");
            out.texts(layout.infections);
            out.text("latent_parents");
            out.texts(layout.latentParents);
            out.text("loci");
            out.texts(layout.loci);
            out.text("samplers");
            out.texts(layout.samplers);
        }

        void write(CborFile& out, const ChainState& state) {
            out.map(9);
            out.text("latent_genotypes");
            out.bytes(state.latentGenotypes);
            out.text("allele_frequencies");
            out.array(state.alleleFrequencies.size());
            for (const auto& frequencies : state.alleleFrequencies) {
                out.numbers(frequencies);
            }
            out.text("ordering");
            out.array(state.ordering.size());
            for (const auto at : state.ordering) {
                out.number(std::uint64_t{at});
            }
            out.text("infection_durations");
            out.numbers(state.infectionDurations);
            out.text("expected_false_positives");
            out.numbers(state.expectedFalsePositives);
            out.text("expected_false_negatives");
            out.numbers(state.expectedFalseNegatives);
            out.text("mean_coi");
            out.number(state.meanCoi);
            out.text("mean_strains_tx");
            out.number(state.meanStrainsTransmitted);
            out.text("parent_set_size_prob");
            out.number(state.parentSetSizeProb);
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
            write(out, replica.samplerCounts);
            out.text("adaptation");
            write(out, replica.adaptation);
            out.text("scheduler_iterations");
            out.number(replica.schedulerIterations);
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
            write(out, checkpoint.layout);
            out.text("run_rng");
            out.text(checkpoint.runRng);
            out.text("replicas");
            out.array(checkpoint.replicas.size());
            for (const auto& replica : checkpoint.replicas) {
                write(out, replica);
            }
            out.text("exchange");
            out.value(toJson(checkpoint.exchange));
        }

        Checkpoint fromJson(const json& j) {
            Checkpoint checkpoint{.identity     = identityFromJson(j.at("identity")),
                                  .position     = positionFromJson(j.at("run_position")),
                                  .parentSetIds = parentSetIdsFromJson(j.at("parent_set_ids")),
                                  .layout       = {.infections    = j.at("layout").at("infections"),
                                                   .latentParents = j.at("layout").at("latent_parents"),
                                                   .loci          = j.at("layout").at("loci"),
                                                   .samplers      = j.at("layout").at("samplers")}};
            checkpoint.runRng = j.at("run_rng").get<std::string>();
            for (const auto& replica : j.at("replicas")) {
                checkpoint.replicas.push_back(replicaFromJson(replica));
            }
            checkpoint.exchange = exchangeFromJson(j.at("exchange"));
            return checkpoint;
        }

    }// namespace

    Checkpoint capture(const Ensemble& ensemble, const Engine& runRng, Identity identity, RunPosition position, ParentSetIdDictionary parentSetIds) {
        Checkpoint checkpoint{.identity = std::move(identity), .position = position, .parentSetIds = std::move(parentSetIds), .runRng = engineState(runRng), .exchange = ensemble.exchangeState()};
        checkpoint.layout          = stateLayout(*ensemble.chains.front().state);
        checkpoint.layout.samplers = ensemble.chains.front().sampler->scheduler_.samplerIds();
        const auto& ids = checkpoint.layout.samplers;
        for (const auto& chain : ensemble.chains) {
            const auto& scheduler = chain.sampler->scheduler_;
            if (scheduler.samplerIds() != ids or infectionIds(chain.state->infections) != checkpoint.layout.infections) {
                throw std::logic_error("chains of one ensemble have different samplers or Infections");
            }
            Replica replica{.state = captureState(*chain.state), .rng = engineState(*chain.r), .schedulerIterations = scheduler.iterations()};
            replica.samplerCounts.reserve(ids.size());
            for (std::uint32_t ii = 0; ii < ids.size(); ++ii) {
                const auto& sampler = scheduler.sampler(ids[ii]);
                replica.samplerCounts.push_back(sampler.counts());
                if (auto adaptation = sampler.adaptationState(); !adaptation.scales.empty() or adaptation.updates != 0) {
                    replica.adaptation.emplace(ii, std::move(adaptation));
                }
            }
            checkpoint.replicas.push_back(std::move(replica));
        }
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
        const auto numChains = static_cast<int>(checkpoint.replicas.size());
        auto ensemble = std::make_unique<Ensemble>(numChains, 0.5, runRng, numCores, [&](const size_t replica, std::shared_ptr<Engine> chainRng) {
            auto state = std::make_shared<State>(input, symptomaticIDPDist, asymptomaticIDPDist, chainRng, nullModel);
            restoreState(checkpoint.layout, checkpoint.replicas[replica].state, *state);
            return state;
        });

        // The burn-in iteration count isn't saved: in burn-in it is the position's iteration. Once burn-in has ended
        // it is never read.
        auto exchange             = checkpoint.exchange;
        exchange.burnInIterations = static_cast<int>(checkpoint.position.iteration);
        // Once burn-in has ended, barriers still pending can't be told apart by pair: the ladder works out each pair's
        // attempts from burn-in iterations that run back to back, and Sampling iterations have alternated the pairs
        // since. A Run started from here begins its own.
        if (checkpoint.position.phase == Phase::Sampling) {
            std::ranges::fill(exchange.swapBarriers, 0.0);
            exchange.barrierIterations = 0;
        }
        ensemble->restoreExchangeState(exchange);
        std::map<std::string, std::uint32_t> samplerAt;// position of each sampler id in the checkpoint's Layout
        for (std::uint32_t ii = 0; ii < checkpoint.layout.samplers.size(); ++ii) {
            samplerAt.emplace(checkpoint.layout.samplers[ii], ii);
        }
        for (size_t replica = 0; replica < ensemble->chains.size(); ++replica) {
            const auto& saved = checkpoint.replicas[replica];
            auto& chain = ensemble->chains[replica];
            setEngineState(*chain.r, saved.rng);
            auto& scheduler = chain.sampler->scheduler_;
            // The scheduler counts in int.
            if (saved.schedulerIterations < 0 or saved.schedulerIterations > std::numeric_limits<int>::max()) {
                refuse(fmt::format("checkpoint of Run {} has replica {} at scheduler iteration {}, out of range", checkpoint.identity.runId, replica, saved.schedulerIterations));
            }
            scheduler.restoreIterations(static_cast<int>(saved.schedulerIterations));
            if (saved.samplerCounts.size() != samplerAt.size()) {
                refuse(fmt::format("checkpoint of Run {} has counts for {} samplers in replica {}, not {}", checkpoint.identity.runId, saved.samplerCounts.size(), replica, samplerAt.size()));
            }
            for (const auto& id : scheduler.samplerIds()) {
                const auto at = savedFor(samplerAt, id, "sampler");
                auto& sampler = scheduler.sampler(id);
                sampler.restoreCounts(saved.samplerCounts[at]);
                if (const auto adaptation = saved.adaptation.find(at); adaptation != saved.adaptation.end()) {
                    sampler.restoreAdaptationState(adaptation->second);
                }
            }
        }
        setEngineState(*runRng, checkpoint.runRng);
        return ensemble;
    }
}// namespace transmission_nets::impl::Model::checkpoint
