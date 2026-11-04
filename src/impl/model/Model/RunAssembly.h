//
// Run assembly: turns a resolved Run recipe, the loaded inputs and the checkpoint if there is one into the pieces a Run
// runs as: the top-level engine, the Ensemble, the Run recorder and the Run driver, with where the Run starts. The one
// place that builds a Run, whether new, resumed or started from a checkpoint; main and the tests both call it. It owns
// the seed draw, the StartFrom reseed and burn-in reset, and the Run manifest's identity. It reads no files: the store
// is whatever sink it is given.
//

#ifndef TRANSMISSION_NETWORKS_APP_RUNASSEMBLY_H
#define TRANSMISSION_NETWORKS_APP_RUNASSEMBLY_H

#include "Checkpoint.h"
#include "Ensemble.h"
#include "Progress.h"
#include "RunDriver.h"
#include "RunRecipe.h"
#include "RunRecorder.h"
#include "Study.h"

#include "core/io/StoreSink.h"

#include <chrono>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace transmission_nets::impl::Model {
    // What only the caller knows about a new Run, for its Run manifest. A resumed Run keeps its stored manifest and
    // takes nothing from this.
    struct Provenance {
        std::string runId;
        std::vector<std::string> cliArguments{};
        std::string gitCommit;
    };

    // 128 random bits as hex, drawn apart from the Run's seeded engines so that it changes nothing they draw.
    std::string newRunId();

    // A Run assembled from its parts, ready to drive. The driver holds references into the other members, so the
    // struct may move but its members may not be replaced.
    struct AssembledRun {
        std::shared_ptr<checkpoint::Engine> rng;
        std::unique_ptr<Ensemble> ensemble;
        std::unique_ptr<RunRecorder> recorder;
        std::unique_ptr<RunDriver> driver;
        RunStart start;
        RunPlan plan;
        RunManifest manifest;// as written to the store: a new Run's with its seed, Run id and start time filled in

        // Runs the Run from where it starts to the end of its plan, or until `stop` asks it to stop. Call once.
        RunOutcome drive(const RunDriver::StopQuery& stop) {
            return driver->drive(start, plan, stop);
        }
    };

    // Assembles the Run `recipe` describes, of `study`. `saved` is the checkpoint the recipe resumes or starts from,
    // none for a new Run. `sink` takes the Run's store; the Run's Checkpoints go to the recipe's store directory. `clock` gives the
    // manifest's start time and the Run recorder's and Run driver's wall times.
    AssembledRun assemble(const Recipe& recipe, const Study& study, const std::optional<checkpoint::Checkpoint>& saved,
                          std::unique_ptr<core::io::StoreSink> sink, Progress& progress, const Provenance& provenance,
                          RunRecorder::Clock clock = std::chrono::system_clock::now);
}// namespace transmission_nets::impl::Model

#endif//TRANSMISSION_NETWORKS_APP_RUNASSEMBLY_H
