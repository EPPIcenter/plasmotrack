//
// The Run recorder writing a store. See RunRecorder.h.
//

#include "RunRecorder.h"

#include "core/utils/timers.h"

#include <fmt/core.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>

namespace transmission_nets::impl::Model {
    using core::io::Table;
    using core::utils::timers::wallTime;

    namespace {
        std::string phaseName(const Phase phase) {
            return phase == Phase::BurnIn ? "burn_in" : "sampling";
        }

        // Appends one row; the values are given in column order and must match the columns' types exactly.
        template<typename... Values>
        void addRow(Table& table, const Values&... values) {
            std::size_t k = 0;
            (std::get<std::vector<Values>>(table.columns.at(k++).second).push_back(values), ...);
        }

        using I = std::vector<std::int64_t>;
        using D = std::vector<double>;
        using S = std::vector<std::string>;

        // The leading columns of every Trace table.
        std::vector<std::pair<std::string, core::io::Column>> traceColumns(std::vector<std::pair<std::string, core::io::Column>> columns) {
            columns.insert(columns.begin(), {{"iteration", I{}}, {"phase", S{}}, {"window", I{}}});
            return columns;
        }

        double seconds(const std::chrono::system_clock::duration d) {
            return std::chrono::duration<double>(d).count();
        }
    }// namespace

    // Rows of every table for the in-flight Window not yet handed to the sink. The run_scalars and infections tables
    // have one column per scalar and per-Infection Quantity, in the State's enumeration order.
    struct RunRecorder::WindowRows {
        explicit WindowRows(const State& structure) {
            std::vector<std::pair<std::string, core::io::Column>> scalars;
            for (const auto& name : structure.scalarNames()) {
                scalars.emplace_back(name, D{});
            }
            for (const auto& term : {"prior", "likelihood", "observation", "transmission", "source"}) {
                scalars.emplace_back(term, D{});
            }
            tables.emplace("run_scalars", Table{traceColumns(std::move(scalars))});
            std::vector<std::pair<std::string, core::io::Column>> infections{{"infection", S{}}};
            for (const auto& name : structure.perInfectionNames()) {
                infections.emplace_back(name, D{});
            }
            tables.emplace("infections", Table{traceColumns(std::move(infections))});
        }

        std::map<std::string, Table> tables{
                {"ordering", {traceColumns({{"infection", S{}}, {"rank", I{}}})}},
                {"allele_frequencies", {traceColumns({{"locus", S{}}, {"allele", I{}}, {"frequency", D{}}})}},
                {"latent_genotypes", {traceColumns({{"infection", S{}}, {"locus", S{}}, {"alleles", S{}}})}},
                {"latent_parent_genotypes", {traceColumns({{"latent_parent", S{}}, {"locus", S{}}, {"alleles", S{}}})}},
                {"parent_sets", {traceColumns({{"infection", S{}}, {"set_id", I{}}, {"prob", D{}}})}},
                // Blocks: counts over the Window's Sampling iterations, added when it closes.
                {"acceptance", {{{"window", I{}}, {"phase", S{}}, {"sampler", S{}}, {"accepted", I{}}, {"proposed", I{}}}}},
                {"swaps", {{{"window", I{}}, {"phase", S{}}, {"temperature_pair", I{}}, {"accepted", I{}}, {"attempted", I{}}}}},
                {"events", {traceColumns({{"wall_time", S{}}, {"type", S{}}, {"payload", S{}}})}},
                {"chains", {traceColumns({{"temperature_slot", I{}}, {"replica", I{}}, {"temperature", D{}}, {"prior", D{}}, {"likelihood", D{}}})}},
                // A lookup table, not a Trace: the members of each Parent set new in this Window.
                {"parent_set_members", {{{"infection", S{}}, {"set_id", I{}}, {"member", S{}}}}},
        };
    };

    namespace {
        nlohmann::json manifestJson(const RunManifest& m, const RecorderOptions& options) {
            nlohmann::json j{
                    {"run_id", m.runId},
                    {"started_at", m.startedAt},
                    {"seed", m.seed},
                    {"cli_arguments", m.cliArguments},
                    {"engine_version", m.engineVersion},
                    {"git_commit", m.gitCommit},
                    {"input_sha256", m.inputSha256},
                    {"chains", m.chains},
                    {"burn_in", m.burnIn},
                    {"samples", m.samples},
                    {"null_model", m.nullModel},
                    {"thin", options.thin},
                    {"window", options.window},
                    {"parent_set_floor", options.parentSetFloor},
            };
            if (m.gradient) {
                j["gradient"] = *m.gradient;
            }
            if (m.startedFrom) {
                j["started_from"] = *m.startedFrom;
            }
            return j;
        }
    }// namespace

    StoredRun storedRunFromJson(const nlohmann::json& j) {
        RunManifest m{.runId         = j.at("run_id").get<std::string>(),
                      .startedAt     = j.at("started_at").get<std::string>(),
                      .seed          = j.at("seed").get<long>(),
                      .cliArguments  = j.at("cli_arguments").get<std::vector<std::string>>(),
                      .engineVersion = j.at("engine_version").get<std::string>(),
                      .gitCommit     = j.at("git_commit").get<std::string>(),
                      .inputSha256   = j.at("input_sha256").get<std::map<std::string, std::string>>(),
                      .chains        = j.at("chains").get<int>(),
                      .burnIn        = j.at("burn_in").get<std::int64_t>(),
                      .samples       = j.at("samples").get<std::int64_t>(),
                      .nullModel     = j.at("null_model").get<bool>()};
        if (j.contains("gradient")) {
            m.gradient = j.at("gradient").get<double>();
        }
        if (j.contains("started_from")) {
            m.startedFrom = j.at("started_from");
        }
        return {.manifest = std::move(m),
                .options  = {.thin = j.at("thin").get<std::int64_t>(), .window = j.at("window").get<std::int64_t>(), .parentSetFloor = j.at("parent_set_floor").get<double>()}};
    }

    std::optional<StoredRun> readStoredRun(const std::filesystem::path& store) {
        const auto file = store / "manifest.json";
        if (!std::filesystem::exists(file)) {
            return std::nullopt;
        }
        try {
            std::ifstream in(file);
            return storedRunFromJson(nlohmann::json::parse(in));
        } catch (const nlohmann::json::exception& e) {
            throw Refused({.kind = RefusalKind::UnreadableRunManifest, .message = fmt::format("{} is not a readable Run manifest: {}", file.string(), e.what())});
        }
    }

    RunRecorder::RunRecorder(const State& structure, std::vector<std::string> samplerIds, std::unique_ptr<core::io::StoreSink> sink, const RecorderOptions options, const RunManifest& manifest, Clock clock)
        : samplerIds_(std::move(samplerIds)), chains_(manifest.chains), sink_(std::move(sink)), options_(options), clock_(std::move(clock)),
          burnIn_(manifest.burnIn), samples_(manifest.samples), rows_(std::make_unique<WindowRows>(structure)) {
        for (const auto& infection : structure.infections) {
            infectionIds_.push_back(infection->id());
        }
        for (const auto& latentParent : structure.latentParents) {
            latentParentIds_.push_back(latentParent->id());
        }
        for (const auto& [label, _] : structure.loci) {
            locusLabels_.push_back(label);
        }

        sink_->writeManifest(manifestJson(manifest, options_));
        writeStructure(structure);
    }

    // Lookup tables a reader needs to decode the Traces and Blocks: Allele positions per Locus, each Infection's Latent
    // parent, and each sampler's type and target.
    void RunRecorder::writeStructure(const State& structure) {
        std::vector<std::string> locus;
        std::vector<std::int64_t> allele;
        for (const auto& [label, l] : structure.loci) {
            for (unsigned int a = 0; a < l->totalAlleles(); ++a) {
                locus.push_back(label);
                allele.push_back(a);
            }
        }
        sink_->writeStructure("loci", {.columns = {{"locus", std::move(locus)}, {"allele", std::move(allele)}}});

        std::vector<bool> symptomatic;
        for (const auto& infection : structure.infections) {
            symptomatic.push_back(infection->isSymptomatic());
        }
        sink_->writeStructure("infections", {.columns = {{"infection", infectionIds_}, {"latent_parent", latentParentIds_}, {"symptomatic", std::move(symptomatic)}}});

        std::vector<std::string> type;
        std::vector<std::string> target;
        for (const auto& id : samplerIds_) {
            const auto slash = id.find('/');
            type.push_back(id.substr(0, slash));
            target.push_back(slash == std::string::npos ? "" : id.substr(slash + 1));
        }
        sink_->writeStructure("samplers", {.columns = {{"sampler", samplerIds_}, {"type", std::move(type)}, {"target", std::move(target)}}});
    }

    RunRecorder::RunRecorder(RunRecorder&&) noexcept = default;
    RunRecorder::~RunRecorder()                        = default;

    void RunRecorder::onWindowCommitted(WindowCommitted committed) {
        windowCommitted_ = std::move(committed);
    }

    // Between Windows: `window_` is one past the last committed, and `phase_` that of the next Sampling iteration.
    RunPosition RunRecorder::position() const {
        const bool burningIn = phase_ == Phase::BurnIn;
        return {.iteration         = iteration_,
                .phase             = phase_,
                .lastClosedWindow  = window_ - 1,
                .burnInRemaining   = burningIn ? burnIn_ - iteration_ : 0,
                .samplingRemaining = burningIn ? samples_ : samples_ - (iteration_ - burnIn_)};
    }

    void RunRecorder::recordStart(const EnsembleSnapshot& ensemble) {
        startedAt_      = clock_();
        windowOpenedAt_ = startedAt_;
        lastRecordedAt_ = startedAt_;
        takeCounts(ensemble);
        blockStart_ = latest_;
        ++windowIterations_;
        raise("run_started", {{"started_at", wallTime(startedAt_)}, {"burn_in", burnIn_}, {"samples", samples_}});
        checkFinite(ensemble);
        writeTraces(ensemble);
    }

    // The `-t` and Window grids count from the Phase's start, which the Run's burn-in length places.
    void RunRecorder::recordResume(const Resumption& resumption, const EnsembleSnapshot& ensemble) {
        const auto& position = resumption.position;
        sink_->discardWindowsAfter(position.lastClosedWindow);
        startedAt_           = clock_();
        windowOpenedAt_      = startedAt_;
        lastRecordedAt_      = startedAt_;
        iteration_           = position.iteration;
        phase_               = position.phase;
        phaseIterations_     = position.phase == Phase::Sampling ? position.iteration - burnIn_ : position.iteration;
        window_              = position.lastClosedWindow + 1;
        samples_             = resumption.samples;
        parentSetIds_        = resumption.parentSetIds;
        takeCounts(ensemble);
        blockStart_ = latest_;
        raise("run_resumed", {{"iteration", position.iteration}, {"window", position.lastClosedWindow}, {"checkpoint", resumption.checkpoint}, {"samples", resumption.samples}});
        checkFinite(ensemble);
    }

    // Trace rows fall on the `-t` grid of the Phase, and Windows on a grid of `window`, both counted from the Phase's
    // start, so every full Window holds window / thin rows per Trace. A complete Window's time ends when its last
    // Sampling iteration was recorded: this one's time belongs to the next Window.
    bool RunRecorder::recordIteration(const EnsembleSnapshot& ensemble) {
        const auto now = clock_();
        if (complete_) {
            commitWindow(lastRecordedAt_);
        }
        lastRecordedAt_ = now;
        ++iteration_;
        ++phaseIterations_;
        ++windowIterations_;
        ++windowSampled_;
        takeCounts(ensemble);
        checkFinite(ensemble);
        const bool written = phaseIterations_ % options_.thin == 0;
        if (written) {
            writeTraces(ensemble);
        }
        complete_ = phaseIterations_ % options_.window == 0;
        return written;
    }

    void RunRecorder::recordTemperatureAdaptation(const TemperatureAdaptation& adaptation) {
        raise("temperatures_adapted", {{"applied", adaptation.applied}, {"temperatures", adaptation.temperatures}, {"barriers", adaptation.barriers}});
    }

    // The proposal variances are by sampler id, one per proposal scale: SALT's are by Allele. A sampler with no
    // proposal scale is left out. The in-flight Window ends here, even if short, so every Window belongs to one Phase.
    void RunRecorder::recordBurnInEnd(const FinalAdaptation& adaptation) {
        if (adaptation.samplers.size() != samplerIds_.size()) {
            throw std::invalid_argument(fmt::format("burn-in ended with adaptation states for {} samplers, not {}", adaptation.samplers.size(), samplerIds_.size()));
        }
        // A Run resumes only from a checkpoint with burn-in left to run, so its first Window always holds some.
        if (windowIterations_ == 0) {
            throw std::logic_error("burn-in ended before any Sampling iteration was recorded in the in-flight Window");
        }
        auto variances = nlohmann::json::object();
        for (std::size_t k = 0; k < samplerIds_.size(); ++k) {
            const auto& scales = adaptation.samplers[k].scales;
            if (scales.empty()) {
                continue;
            }
            auto& v = variances[samplerIds_[k]] = nlohmann::json::array();
            for (const auto& scale : scales) {
                v.push_back(scale.variance);
            }
        }
        raise("burn_in_ended", {{"proposal_variances", std::move(variances)}, {"temperatures", adaptation.temperatures}});
        burnInEnded_ = true;
        complete_    = true;
    }

    void RunRecorder::recordInterrupt(const std::string& signal, const Interrupt which) {
        if (ended_) {
            return;
        }
        interrupted_ = true;
        raise("interrupted", {{"signal", signal}, {"second", which == Interrupt::Second}});
        if (which == Interrupt::Second) {
            forced_ = true;
            recordEnd();
        }
    }

    void RunRecorder::recordEnd() {
        if (ended_) {
            return;
        }
        ended_         = true;
        const auto now = clock_();
        raise("run_ended", {{"reason", interrupted_ ? "interrupted" : "completed"}, {"seconds", seconds(now - startedAt_)}});
        if (windowIterations_ > 0) {
            commitWindow(now);
        }
    }

    // A Run event at the last Sampling iteration recorded, in the in-flight Window.
    void RunRecorder::raise(const std::string& type, nlohmann::json payload) {
        addRow(rows_->tables.at("events"), iteration_, phaseName(phase_), window_, wallTime(clock_()), type, payload.dump());
    }

    // One row per Trace table for the Hot chain as it is now. Rows read the Hot chain by position (Infection, Latent
    // parent) or label (Locus), never by object, because each chain has its own objects.
    void RunRecorder::writeTraces(const EnsembleSnapshot& ensemble) {
        State& s                  = ensemble.hot.state;
        Model& m                  = ensemble.hot.model;
        auto& t                   = rows_->tables;
        const std::string phase   = phaseName(phase_);
        const std::int64_t it     = iteration_;
        const std::int64_t window = window_;

        double observation = 0;
        for (const auto& process : m.observationProcessLikelihoodList) {
            observation += process->value();
        }
        double transmission = 0;
        for (const auto& process : m.transmissionProcessList) {
            transmission += process->value();
        }
        // Kept for the table's layout: latent parents are summed out of the Transmission process, so there is no
        // separate Source term.
        const double source = 0;

        // The Parent-set distribution of each Infection, normalised over every Parent set, then floored. A Parent set with
        // probability exactly 0 (impossible, or underflowed) carries nothing, so it is dropped even with a floor of 0
        // (#67). An Infection's rows go out in ascending set id: the scoring order shifts as candidate parents come and go, while sorted
        // rows repeat from one iteration to the next and compress (#55).
        std::vector<std::pair<std::int64_t, double>> infectionSets;
        for (std::size_t i = 0; i < m.transmissionProcessList.size(); ++i) {
            infectionSets.clear();
            const auto dist = m.transmissionProcessList[i]->calcParentSetDist();
            for (const auto& [llik, ps] : dist.parentSetLliks) {
                const double prob = std::exp(llik - dist.totalLlik);
                if (prob == 0 or prob < options_.parentSetFloor) {
                    continue;
                }
                std::vector<std::string> members;
                for (const auto& member : ps) {
                    members.push_back(member->id());
                }
                std::ranges::sort(members);
                auto& ids                  = parentSetIds_[infectionIds_[i]];
                const auto [found, isNew] = ids.ids.try_emplace(members, ids.nextId);
                const auto setId          = static_cast<std::int64_t>(found->second);
                if (isNew) {
                    ++ids.nextId;
                    for (const auto& member : members) {
                        addRow(t.at("parent_set_members"), infectionIds_[i], setId, member);
                    }
                }
                infectionSets.emplace_back(setId, prob);
            }
            std::ranges::sort(infectionSets);
            for (const auto& [setId, prob] : infectionSets) {
                addRow(t.at("parent_sets"), it, phase, window, infectionIds_[i], setId, prob);
            }
        }

        // Every Quantity's rows. Scalars and per-Infection values are gathered first: a run_scalars row carries every
        // scalar beside the Model's terms, and an infections row carries every per-Infection Quantity.
        struct Tracing {
            RunRecorder& recorder;
            std::map<std::string, Table>& tables;
            std::int64_t it;
            std::string phase;
            std::int64_t window;
            std::vector<double> scalars{};
            std::vector<std::vector<double>> perInfectionValues{};// by Infection, then by Quantity in enumeration order

            void scalar(const char*, const std::shared_ptr<core::parameters::Parameter<double>>& p) {
                scalars.push_back(p->value());
            }
            void perInfection(const char*, const size_t index, const std::shared_ptr<core::parameters::Parameter<double>>& p) {
                if (perInfectionValues.size() <= index) {
                    perInfectionValues.resize(index + 1);
                }
                perInfectionValues[index].push_back(p->value());
            }
            void alleleFrequencies(const std::shared_ptr<LocusImpl>& locus, const std::shared_ptr<core::parameters::Parameter<AlleleFrequencyImpl>>& p) {
                const auto& frequencies = p->value();
                for (unsigned int a = 0; a < frequencies.totalElements(); ++a) {
                    addRow(tables.at("allele_frequencies"), it, phase, window, locus->label, static_cast<std::int64_t>(a), frequencies.frequencies(a));
                }
            }
            void latentGenotype(const size_t index, const std::shared_ptr<LocusImpl>& locus, const std::shared_ptr<core::parameters::Parameter<GeneticsImpl>>& p) {
                addRow(tables.at("latent_genotypes"), it, phase, window, recorder.infectionIds_[index], locus->label, p->value().serialize());
            }
            void ordering(const std::shared_ptr<OrderingImpl>& ordering) {
                const auto order = ordering->value();
                for (std::size_t rank = 0; rank < order.size(); ++rank) {
                    addRow(tables.at("ordering"), it, phase, window, order[rank]->id(), static_cast<std::int64_t>(rank));
                }
            }
        };
        Tracing tracing{.recorder = *this, .tables = t, .it = it, .phase = phase, .window = window};
        s.quantities(tracing);
        // A row given column by column, for the tables whose columns the State's enumeration decides.
        const auto push = [](Table& table, std::size_t& k, const auto& value) {
            std::get<std::vector<std::decay_t<decltype(value)>>>(table.columns.at(k++).second).push_back(value);
        };
        {
            auto& table   = t.at("run_scalars");
            std::size_t k = 0;
            push(table, k, it);
            push(table, k, phase);
            push(table, k, window);
            for (const double value : tracing.scalars) {
                push(table, k, value);
            }
            for (const double term : {m.getPrior(), m.getLikelihood(), observation, transmission, source}) {
                push(table, k, term);
            }
        }
        for (std::size_t i = 0; i < infectionIds_.size(); ++i) {
            auto& table   = t.at("infections");
            std::size_t k = 0;
            push(table, k, it);
            push(table, k, phase);
            push(table, k, window);
            push(table, k, infectionIds_[i]);
            for (const double value : tracing.perInfectionValues.at(i)) {
                push(table, k, value);
            }
        }
        // latent_parent_genotypes stays empty: latent parents are summed out and have no genotypes to record.

        for (std::size_t slot = 0; slot < ensemble.slots.size(); ++slot) {
            const auto& [replica, temperature, prior, likelihood] = ensemble.slots[slot];
            addRow(t.at("chains"), it, phase, window, static_cast<std::int64_t>(slot), static_cast<std::int64_t>(replica), temperature, prior, likelihood);
        }

        handOverRows();
        sink_->allowRowGroupEnd();
    }

    // Hands the rows of every table so far to the sink and empties them, keeping their room for the next Trace row.
    void RunRecorder::handOverRows() {
        for (auto& [name, table] : rows_->tables) {
            if (table.rows() > 0) {
                sink_->appendRows(name, table);
                for (auto& [_, column] : table.columns) {
                    std::visit([](auto& values) { values.clear(); }, column);
                }
            }
        }
    }

    // A replica's counts since the last call are the Hot chain's only if it sat in temperature slot 0 while that
    // Sampling iteration ran. Swaps come last in a Sampling iteration, so that is the replica in slot 0 at the last
    // call, not at this one. Replicas move between slots, so one replica's counts over a Window would mix in Sampling
    // iterations it ran at other temperatures.
    void RunRecorder::takeCounts(const EnsembleSnapshot& ensemble) {
        const bool fits = ensemble.slots.size() == chains_ and ensemble.swaps.size() == chains_ - 1 and ensemble.acceptance.size() == chains_ and
                          std::ranges::all_of(ensemble.acceptance, [&](const auto& counts) { return counts.size() == samplerIds_.size(); });
        if (!fits) {
            throw std::invalid_argument(fmt::format("an ensemble of {} chains and {} samplers needs one temperature slot per chain, one count per sampler per replica and one per neighbouring pair of slots; got {} slots, {} replicas' counts and {} pairs",
                                                    chains_, samplerIds_.size(), ensemble.slots.size(), ensemble.acceptance.size(), ensemble.swaps.size()));
        }
        if (replicaCounts_.empty()) {
            latest_.acceptance.assign(samplerIds_.size(), {});
        } else {
            const auto& before = replicaCounts_[hotReplica_];
            const auto& now    = ensemble.acceptance[hotReplica_];
            for (std::size_t k = 0; k < samplerIds_.size(); ++k) {
                latest_.acceptance[k].accepted += now[k].accepted - before[k].accepted;
                latest_.acceptance[k].proposed += now[k].proposed - before[k].proposed;
            }
        }
        replicaCounts_ = ensemble.acceptance;
        hotReplica_    = ensemble.slots.front().replica;
        latest_.swaps  = ensemble.swaps;
        hot_           = {.prior = ensemble.slots.front().prior, .likelihood = ensemble.slots.front().likelihood};
    }

    // Checks every temperature slot's prior and likelihood at every call, thinned or not. A value is raised when it
    // becomes NaN or infinite, not again while it stays so.
    void RunRecorder::checkFinite(const EnsembleSnapshot& ensemble) {
        nonFinite_.resize(ensemble.slots.size());
        auto check = [&](const std::size_t slot, const char* term, const double value, bool& wasNonFinite) {
            const bool nonFinite = !std::isfinite(value);
            if (nonFinite and !wasNonFinite) {
                raise("non_finite_value", {{"temperature_slot", slot}, {"term", term}, {"value", std::isnan(value) ? "nan" : value > 0 ? "inf" : "-inf"}});
            }
            wasNonFinite = nonFinite;
        };
        for (std::size_t slot = 0; slot < ensemble.slots.size(); ++slot) {
            check(slot, "prior", ensemble.slots[slot].prior, nonFinite_[slot].first);
            check(slot, "likelihood", ensemble.slots[slot].likelihood, nonFinite_[slot].second);
        }
    }

    // One row per sampler and per neighbouring pair of temperature slots: what the counts gained since the last Block.
    // Returns the swap counts written.
    std::vector<SwapCounts> RunRecorder::writeBlocks() {
        auto& t                 = rows_->tables;
        const std::string phase = phaseName(phase_);
        auto gained             = [](const std::uint64_t now, const std::uint64_t before) { return static_cast<std::int64_t>(now - before); };
        for (std::size_t k = 0; k < samplerIds_.size(); ++k) {
            const auto& [accepted, proposed] = latest_.acceptance[k];
            const auto& before               = blockStart_.acceptance[k];
            addRow(t.at("acceptance"), window_, phase, samplerIds_[k], gained(accepted, before.accepted), gained(proposed, before.proposed));
        }
        std::vector<SwapCounts> swaps;
        for (std::size_t pair = 0; pair < latest_.swaps.size(); ++pair) {
            const auto& [accepted, attempted] = latest_.swaps[pair];
            const auto& before                = blockStart_.swaps[pair];
            swaps.push_back({accepted - before.accepted, attempted - before.attempted});
            addRow(t.at("swaps"), window_, phase, static_cast<std::int64_t>(pair), gained(accepted, before.accepted), gained(attempted, before.attempted));
        }
        blockStart_ = latest_;
        return swaps;
    }

    // Hands the rest of the in-flight Window's rows to the sink and commits it, then reports where the Run
    // stands. The Window covers at least one Sampling iteration; its time ends at `closedAt`. If burn-in end is
    // recorded in it, the sampling Phase starts with the next Window.
    void RunRecorder::commitWindow(const std::chrono::system_clock::time_point closedAt) {
        WindowSummary summary{.window             = window_,
                              .phase              = phase_,
                              .samplingIterations = windowSampled_,
                              .seconds            = seconds(closedAt - windowOpenedAt_),
                              .hotPrior           = hot_.prior,
                              .hotLikelihood      = hot_.likelihood,
                              .swaps              = std::vector<SwapCounts>(chains_ - 1)};
        if (windowSampled_ > 0) {
            summary.swaps = writeBlocks();
        }
        const double spent = summary.seconds;
        raise("window_closed", {{"sampling_iterations", windowSampled_},
                                {"seconds", spent},
                                {"iterations_per_second", spent > 0 ? nlohmann::json(static_cast<double>(windowSampled_) / spent) : nlohmann::json(nullptr)}});
        windowOpenedAt_ = closedAt;
        handOverRows();
        sink_->commitWindow({.window = window_, .firstIteration = iteration_ - windowIterations_ + 1, .lastIteration = iteration_, .phase = phaseName(phase_)});
        ++window_;
        windowIterations_ = 0;
        windowSampled_    = 0;
        complete_         = false;

        const bool burnInEnded = burnInEnded_;
        if (burnInEnded) {
            phase_           = Phase::Sampling;
            phaseIterations_ = 0;
            burnInEnded_     = false;
        }
        if (windowCommitted_ and !forced_) {
            windowCommitted_({.position = position(), .parentSetIds = parentSetIds_, .burnInEnded = burnInEnded, .summary = summary});
        }
    }
}// namespace transmission_nets::impl::Model
