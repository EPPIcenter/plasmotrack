//
// The Run recipe. See RunRecipe.h.
//

#include "RunRecipe.h"

#include "core/version.h"

#include <fmt/core.h>

namespace transmission_nets::impl::Model {

    namespace {
        std::optional<Refusal> modeRefusal(const GivenOptions& given, const std::optional<StoredRun>& stored) {
            if (given.resume and given.startFrom) {
                return refuse(RefusalKind::ResumeWithStartFrom, "--resume and --start-from can't be given together.");
            }
            if (given.resume and given.outputDir) {
                return refuse(RefusalKind::ResumeWithOutputDir, "--resume continues the Run in its own store; give no output directory.");
            }
            if (given.resume and !stored) {
                return refuse(RefusalKind::NothingToResume, fmt::format("{} holds no Run to resume.", given.resume->string()));
            }
            if (!given.resume and !given.outputDir) {
                return refuse(RefusalKind::NoOutputDir, "the option '--output-dir' is required but missing");
            }
            if (!given.resume and stored) {
                return refuse(RefusalKind::OutputHoldsRun, fmt::format("Output directory {} already holds a Run.", given.outputDir->string()));
            }
            return std::nullopt;
        }

        std::optional<Refusal> rangeRefusal(const RecorderOptions& options) {
            if (options.thin < 1) {
                return refuse(RefusalKind::ThinOutOfRange, "Thin must be at least 1.");
            }
            if (options.window < 1 or options.window % options.thin != 0) {
                return refuse(RefusalKind::WindowOutOfRange, fmt::format("--window must be a positive multiple of -t (got --window {}, -t {}).", options.window, options.thin));
            }
            if (!(options.parentSetFloor >= 0 and options.parentSetFloor < 1)) {
                return refuse(RefusalKind::ParentSetFloorOutOfRange, fmt::format("--parent-set-floor must be at least 0 and below 1 (got {}).", options.parentSetFloor));
            }
            return std::nullopt;
        }

        // Both uses of a checkpoint need the same input and chain count; resuming also needs the same engine version.
        std::optional<Refusal> compatibilityRefusal(const checkpoint::Checkpoint& saved, const Mode mode, const std::string& inputHash, const std::size_t chains) {
            const auto& id = saved.identity;
            if (id.inputHash != inputHash) {
                return refuse(RefusalKind::InputChanged, fmt::format("checkpoint of Run {} was made from input {}, not {}", id.runId, id.inputHash, inputHash));
            }
            if (saved.ensemble.replicas.size() != chains) {
                return refuse(RefusalKind::ChainCountChanged, fmt::format("checkpoint of Run {} has {} chains, not {}", id.runId, saved.ensemble.replicas.size(), chains));
            }
            if (mode == Mode::Resume and id.engineVersion != core::version::getVersionString()) {
                return refuse(RefusalKind::EngineVersionChanged, fmt::format("checkpoint of Run {} was made by engine version {}; resuming needs the same version, not {}", id.runId, id.engineVersion, core::version::getVersionString()));
            }
            return std::nullopt;
        }

        // Sampling iterations the Run planned when `position` was taken.
        std::int64_t plannedSamples(const RunPosition& position, const std::int64_t burnIn) {
            return position.phase == Phase::Sampling ? position.iteration - burnIn + position.samplingRemaining : position.samplingRemaining;
        }

        // The Run in the store, as it was started: an option given again must be what the Run started with, and only
        // -s may change, and only grow.
        RecipeOrRefusal resumedRun(const GivenOptions& given, const StoredRun& stored, const checkpoint::Checkpoint& saved) {
            const auto& manifest = stored.manifest;
            std::optional<Refusal> refusal;
            const auto requireSame = [&](const char* flag, const auto& given, const auto& started) {
                if (!refusal and given and *given != started) {
                    refusal = refuse(RefusalKind::RecipeChanged,
                                     fmt::format("the Run was started with {} {}; resuming it needs the same, not {}.", flag, nlohmann::json(started).dump(), *given));
                }
            };
            requireSame("-n", given.chains, manifest.chains);
            requireSame("-t", given.thin, stored.options.thin);
            requireSame("--window", given.window, stored.options.window);
            requireSame("-b", given.burnIn, manifest.burnIn);
            if (manifest.gradient) {
                requireSame("-g", given.gradient, *manifest.gradient);
            } else if (given.gradient and !refusal) {
                refusal = refuse(RefusalKind::RecipeChanged, "the Run took its temperature ladder from the checkpoint it started from; resuming it takes no -g.");
            }
            requireSame("--parent-set-floor", given.parentSetFloor, stored.options.parentSetFloor);
            requireSame("--seed", given.seed, manifest.seed);
            if (given.nullModel and !manifest.nullModel and !refusal) {
                refusal = refuse(RefusalKind::RecipeChanged, "the Run was started without --null-model; resuming it needs the same.");
            }
            if (refusal) {
                return *refusal;
            }
            for (const auto* idp : {input_names::SYMPTOMATIC_IDP, input_names::ASYMPTOMATIC_IDP}) {
                const auto started = manifest.inputSha256.find(idp);
                const auto now     = given.inputSha256.find(idp);
                if (started == manifest.inputSha256.end() or now == given.inputSha256.end() or started->second != now->second) {
                    return refuse(RefusalKind::IdpChanged, fmt::format("the {} given is not the one the Run was started with.", idp));
                }
            }
            if (saved.identity.runId != manifest.runId) {
                return refuse(RefusalKind::RunIdMismatch, fmt::format("{} holds a checkpoint of Run {}, not of the Run there ({}).", given.resume->string(), saved.identity.runId, manifest.runId));
            }
            if (auto refusal = compatibilityRefusal(saved, Mode::Resume, given.inputSha256.at(input_names::INPUT), manifest.chains)) {
                return *refusal;
            }

            const auto planned = plannedSamples(saved.position, manifest.burnIn);
            if (given.samples and *given.samples < planned) {
                return refuse(RefusalKind::SamplesShrunk, fmt::format("the Run plans {} sampling iterations; resuming it may grow -s, not shrink it to {}.", planned, *given.samples));
            }
            // The manifest stays the original recipe; the Run event of the resume records the new -s.
            return Recipe{.mode     = Mode::Resume,
                          .store    = *given.resume,
                          .manifest = manifest,
                          .recorder = stored.options,
                          .plan     = {.burnIn = manifest.burnIn, .samples = given.samples.value_or(planned)},
                          .cores    = given.cores.value_or(defaults::CORES)};
        }

        Recipe newRun(const GivenOptions& given) {
            Recipe recipe{.mode     = Mode::New,
                          .store    = *given.outputDir,
                          .recorder = {.thin           = given.thin.value_or(defaults::THIN),
                                       .window         = given.window.value_or(defaults::WINDOW),
                                       .parentSetFloor = given.parentSetFloor.value_or(defaults::PARENT_SET_FLOOR)},
                          .plan     = {.burnIn = given.burnIn.value_or(defaults::BURN_IN), .samples = given.samples.value_or(defaults::SAMPLES)},
                          .cores    = given.cores.value_or(defaults::CORES)};
            recipe.manifest = {.seed          = given.seed.value_or(defaults::SEED),
                               .engineVersion = core::version::getVersionString(),
                               .inputSha256   = given.inputSha256,
                               .chains        = given.chains.value_or(defaults::CHAINS),
                               .gradient      = given.gradient.value_or(defaults::GRADIENT),
                               .burnIn        = recipe.plan.burnIn,
                               .samples       = recipe.plan.samples,
                               .nullModel     = given.nullModel};
            return recipe;
        }

        // A new Run from the checkpoint's State and ladder. Its IDPs may differ from the source Run's: the manifest
        // records the hashes of those given.
        RecipeOrRefusal startedRun(const GivenOptions& given, const checkpoint::Checkpoint& saved) {
            auto withChains = given;
            withChains.chains = given.chains.value_or(static_cast<int>(saved.ensemble.replicas.size()));
            auto recipe       = newRun(withChains);
            if (auto refusal = compatibilityRefusal(saved, Mode::StartFrom, given.inputSha256.at(input_names::INPUT), recipe.manifest.chains)) {
                return *refusal;
            }
            recipe.mode                 = Mode::StartFrom;
            recipe.manifest.gradient    = std::nullopt;
            recipe.manifest.startedFrom = nlohmann::json{{"run_id", saved.identity.runId}, {"engine_version", saved.identity.engineVersion}, {"temperatures", saved.ensemble.exchange.temperatures}};
            if (given.gradient) {
                recipe.notices.emplace_back("Ignoring -g: the temperature ladder comes from the checkpoint.");
            }
            return recipe;
        }
    }// namespace

    RecipeOrRefusal resolve(const GivenOptions& given, const std::optional<StoredRun>& stored, const std::optional<checkpoint::Checkpoint>& saved) {
        if (auto refusal = modeRefusal(given, stored)) {
            return *refusal;
        }
        if (given.resume and !saved) {
            return refuse(RefusalKind::UnusableCheckpoint, fmt::format("{} holds no checkpoint to resume from.", given.resume->string()));
        }
        if (given.startFrom and !saved) {
            return refuse(RefusalKind::UnusableCheckpoint, fmt::format("there is no checkpoint at {}", given.startFrom->string()));
        }
        auto resolved = given.resume      ? resumedRun(given, *stored, *saved)
                        : given.startFrom ? startedRun(given, *saved)
                                          : newRun(given);
        if (const auto* recipe = std::get_if<Recipe>(&resolved)) {
            if (auto refusal = rangeRefusal(recipe->recorder)) {
                return *refusal;
            }
        }
        return resolved;
    }
}// namespace transmission_nets::impl::Model
