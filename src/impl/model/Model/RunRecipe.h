//
// The Run recipe: what was given on the command line, the stored Run when resuming and the checkpoint when resuming or
// starting from one, resolved into the recipe a Run uses, or one refusal. It holds every rule about which options a Run
// takes and the defaults of those not given; it reads no files.
//

#ifndef TRANSMISSION_NETWORKS_APP_RUNRECIPE_H
#define TRANSMISSION_NETWORKS_APP_RUNRECIPE_H

#include "Checkpoint.h"
#include "Refusal.h"
#include "RunDriver.h"
#include "RunRecorder.h"

#include <cstdint>
#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace transmission_nets::impl::Model {

    // Defaults of options not given. --help shows them.
    namespace defaults {
        inline constexpr std::int64_t BURN_IN       = 5000;
        inline constexpr std::int64_t SAMPLES       = 10000;
        inline constexpr std::int64_t THIN          = 10;
        inline constexpr std::int64_t WINDOW        = 100;
        inline constexpr double PARENT_SET_FLOOR    = 1e-9;
        inline constexpr int CHAINS                 = 1;
        inline constexpr unsigned int CORES         = 1;
        inline constexpr double GRADIENT            = 0;
        inline constexpr long SEED                  = -1;// draw a seed
    }// namespace defaults

    // What was given on the command line; an option not given is empty.
    struct GivenOptions {
        std::optional<std::int64_t> burnIn{};
        std::optional<std::int64_t> samples{};
        std::optional<std::int64_t> thin{};
        std::optional<std::int64_t> window{};
        std::optional<double> parentSetFloor{};
        std::optional<int> chains{};
        std::optional<unsigned int> cores{};
        std::optional<double> gradient{};
        std::optional<long> seed{};
        bool nullModel = false;
        std::optional<std::filesystem::path> outputDir{};
        std::optional<std::filesystem::path> resume{};
        std::optional<std::filesystem::path> startFrom{};
        // The inputs given, by their input_names, as their SHA-256.
        std::map<std::string, std::string> inputSha256{};
    };

    // Whether the Run is a new one, the Run in the store resumed, or a new one started from a checkpoint.
    enum class Mode { New, Resume, StartFrom };

    struct Recipe {
        Mode mode = Mode::New;
        std::filesystem::path store;
        // A resumed Run's is the stored one. A new Run's has no run id, start time, command line or git commit yet, and
        // a seed of -1 if one is to be drawn.
        RunManifest manifest;
        RecorderOptions recorder;
        RunPlan plan;
        unsigned int cores = defaults::CORES;
        std::vector<std::string> notices{};// for the console: what was given but plays no part
    };

    using RecipeOrRefusal = std::variant<Recipe, Refusal>;

    // `stored` is the Run the store holds, if any: the one given to --resume, or the output directory of a new Run.
    // `saved` is the checkpoint resumed or started from.
    RecipeOrRefusal resolve(const GivenOptions& given, const std::optional<StoredRun>& stored, const std::optional<checkpoint::Checkpoint>& saved);
}// namespace transmission_nets::impl::Model

#endif//TRANSMISSION_NETWORKS_APP_RUNRECIPE_H
