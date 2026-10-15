//
// Run assembly. See RunAssembly.h.
//

#include "RunAssembly.h"

#include "core/utils/timers.h"

#include <fmt/core.h>

#include <random>

namespace transmission_nets::impl::Model {
    std::string newRunId() {
        std::random_device device;
        std::string id;
        for (int k = 0; k < 4; ++k) {
            id += fmt::format("{:08x}", device());
        }
        return id;
    }

    AssembledRun assemble(const Recipe& recipe, const Inputs& inputs, const std::optional<checkpoint::Checkpoint>& saved,
                          std::unique_ptr<core::io::StoreSink> sink, Progress& progress, const Provenance& provenance, RunRecorder::Clock clock) {
        AssembledRun run;
        run.plan            = recipe.plan;
        run.manifest        = recipe.manifest;
        auto& manifest      = run.manifest;
        const bool resuming = recipe.mode == Mode::Resume;

        if (manifest.seed == -1) {
            manifest.seed = core::utils::timers::time().time_since_epoch().count();
        }

        run.rng = std::make_shared<checkpoint::Engine>(manifest.seed);
        if (saved) {
            run.ensemble = checkpoint::restore(*saved, run.rng, inputs.input, inputs.symptomaticIdp, inputs.asymptomaticIdp, manifest.nullModel, recipe.cores);
            if (recipe.mode == Mode::StartFrom) {
                // A new Run's burn-in counts from 0, not from the saved position, so its ladder adapts on the same
                // iterations as a fresh Run's.
                auto exchange             = run.ensemble->exchangeState();
                exchange.burnInIterations = 0;
                run.ensemble->restoreExchangeState(exchange);
                // A new Run draws from its own seed, from the saved State on.
                run.rng->seed(static_cast<std::uint32_t>(manifest.seed));
                for (auto& chain : run.ensemble->chains) {
                    chain.r->seed((*run.rng)());
                }
            }
        } else {
            run.ensemble = std::make_unique<Ensemble>(manifest.chains, *manifest.gradient, run.rng, manifest.nullModel, recipe.cores, inputs.input, inputs.symptomaticIdp, inputs.asymptomaticIdp);
        }
        if (!resuming) {
            // A resumed Run keeps its manifest: the original recipe. Its run_resumed event records the new -s.
            manifest.runId        = provenance.runId;
            manifest.startedAt    = core::utils::timers::wallTime(clock());
            manifest.cliArguments = provenance.cliArguments;
            manifest.gitCommit    = provenance.gitCommit;
        }

        run.recorder = std::make_unique<RunRecorder>(*run.ensemble->hotChain().state, run.ensemble->samplerIds(), std::move(sink), recipe.recorder, manifest, clock);
        const checkpoint::Identity identity{.runId = manifest.runId, .inputHash = manifest.inputSha256.at(input_names::INPUT), .engineVersion = manifest.engineVersion};
        run.driver = std::make_unique<RunDriver>(*run.ensemble, *run.recorder, recipe.store, identity, progress, std::move(clock));

        run.start = FreshRun{};
        if (resuming) {
            run.start = ResumedRun{.position = saved->position, .parentSetIds = saved->parentSetIds};
        }
        return run;
    }
}// namespace transmission_nets::impl::Model
