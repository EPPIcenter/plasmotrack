//
// Created by Maxwell Murphy on 11/15/23.
//

#include "core/io/ParquetStoreSink.h"
#include "core/io/parse_json.h"
#include "core/io/utils.h"
#include "core/utils/InterruptWatcher.h"
#include "core/utils/Sha256.h"
#include "core/version.h"
#include "impl/model/Model/Checkpoint.h"
#include "impl/model/Model/Progress.h"
#include "impl/model/Model/RunAssembly.h"
#include "impl/model/Model/RunRecipe.h"

#include <boost/program_options.hpp>

#include <fmt/core.h>

#include <atomic>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <variant>
#include <vector>


using namespace transmission_nets::impl;
namespace checkpoint = transmission_nets::impl::Model::checkpoint;
using namespace transmission_nets::core::io;
using namespace transmission_nets::core::computation;
using namespace transmission_nets::core::utils;
using namespace transmission_nets::core::version;

namespace {
    // Exit codes. An output, Run manifest or checkpoint that can't be written or used exits with OUTPUT_FAILURE (3),
    // naming the path; see OutputFailure.h.
    const size_t SUCCESS                   = 0;
    const size_t ERROR_IN_COMMAND_LINE     = 1;// a bad option, or a Run refused to start or resume
    const size_t ERROR_UNHANDLED_EXCEPTION = 2;

    // The value of option `name` if it was given, not defaulted for --help.
    template<typename T>
    std::optional<T> givenValue(const boost::program_options::variables_map& vm, const char* name) {
        if (!vm.count(name) or vm[name].defaulted()) {
            return std::nullopt;
        }
        return vm[name].as<T>();
    }
}

int main(int argc, char** argv) {
    // Declared before the watcher so they outlive the watcher's thread, which forces the Run's end through the Run
    // driver once there is one.
    std::optional<Model::AssembledRun> run;
    std::atomic<Model::RunDriver*> driving{nullptr};
    // Constructed before any other thread starts, so every thread leaves the interrupt signals to the watcher.
    const InterruptWatcher watcher([&](const std::string& signal) {
        if (auto* const d = driving.load()) {
            d->forceEnd(watcher.firstSignal(), signal);
        }
    });


    try {
        bool null_model;
        std::string input;
        std::string symptomatic_idp_path;
        std::string asymptomatic_idp_path;

        // Defaults are the Run recipe's; boost holds them only to show them in --help.
        namespace po       = boost::program_options;
        namespace defaults = Model::defaults;
        po::options_description desc("Options");
        auto opts = desc.add_options();
        opts("help", "Display help message");
        opts("version", "Display version information");
        opts("burnin,b", po::value<std::int64_t>()->default_value(defaults::BURN_IN), "Number of Sampling iterations in burn-in");
        opts("sample,s", po::value<std::int64_t>()->default_value(defaults::SAMPLES), "Number of Sampling iterations in the sampling Phase");
        opts("thin,t", po::value<std::int64_t>()->default_value(defaults::THIN), "A Trace row is written every -t Sampling iterations");
        opts("window", po::value<std::int64_t>()->default_value(defaults::WINDOW), "Sampling iterations per Window, a multiple of --thin. Each Window is committed to the store whole, so a killed Run loses at most the Window in flight.");
        opts("parent-set-floor", po::value<double>()->default_value(defaults::PARENT_SET_FLOOR), "Parent sets less probable than this are not recorded. A floor of 0 keeps every Parent set with non-zero probability.");
        opts("numchains,n", po::value<int>()->default_value(defaults::CHAINS), "Number of chains to run in replica exchange algorithm.");
        opts("numcores,c", po::value<unsigned int>()->default_value(defaults::CORES), "Number of cores to use in replica exchange algorithm.");
        opts("gradient,g", po::value<double>()->default_value(defaults::GRADIENT), "Lower temperature of gradient to use in replica exchange algorithm");
        opts("seed", po::value<long>()->default_value(defaults::SEED), "Seed used in random number generator. Runs with the same seed, input and options produce the same results, regardless of numcores. A value of -1 indicates generate a random seed.");
        opts("symptomatic-idp", po::value<std::string>(&symptomatic_idp_path)->required(), "file path to Symptomatic IDP");
        opts("asymptomatic-idp", po::value<std::string>(&asymptomatic_idp_path)->required(), "file path to Symptomatic IDP");
        opts("input,i", po::value<std::string>(&input)->required(), "Input file");
        opts("output-dir,o", po::value<std::string>(), "Output directory: the Run's store is written here. Required unless --resume is given.");
        opts("resume", po::value<std::string>(), "Continue the Run whose store is this directory from its latest checkpoint, with the same input. Give no output directory. Recipe options default to the Run's; only -s may change, and only grow.");
        opts("start-from", po::value<std::string>(), "Start a new Run from the State and adaptation in this checkpoint (a Run's checkpoint.cbor or checkpoint-burn-in.cbor), with the same input. The new Run has its own seed and -b burn-in, and numbers iterations from 0.");
        opts("null-model", po::bool_switch(&null_model)->default_value(false), "Run the null model (no genetics)");

        po::positional_options_description p;
        p.add("input", 1);
        p.add("output-dir", 1);

        po::variables_map vm;
        try {
            po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(),
                      vm);// can throw

            /** --help option
                       */
            if (vm.count("help")) {
                std::cout << "Plasmotrack: Transmission Network Inference" << std::endl
                          << desc << std::endl;
                return SUCCESS;
            }

            /** --version option
                       */
            if (vm.count("version")) {
                std::cout << "Plasmotrack version " << getVersionString() << std::endl;
                return SUCCESS;
            }

            po::notify(vm);// throws on error, so do after help in case
                           // there are any problems
        } catch (po::error& e) {
            std::cerr << "ERROR: " << e.what() << std::endl
                      << std::endl;
            std::cerr << desc << std::endl;
            return ERROR_IN_COMMAND_LINE;
        }

        namespace fs = std::filesystem;
        using Model::RefusalKind;

        Model::GivenOptions given{.burnIn         = givenValue<std::int64_t>(vm, "burnin"),
                                  .samples        = givenValue<std::int64_t>(vm, "sample"),
                                  .thin           = givenValue<std::int64_t>(vm, "thin"),
                                  .window         = givenValue<std::int64_t>(vm, "window"),
                                  .parentSetFloor = givenValue<double>(vm, "parent-set-floor"),
                                  .chains         = givenValue<int>(vm, "numchains"),
                                  .cores          = givenValue<unsigned int>(vm, "numcores"),
                                  .gradient       = givenValue<double>(vm, "gradient"),
                                  .seed           = givenValue<long>(vm, "seed"),
                                  .nullModel      = null_model,
                                  .outputDir      = givenValue<std::string>(vm, "output-dir"),
                                  .resume         = givenValue<std::string>(vm, "resume"),
                                  .startFrom      = givenValue<std::string>(vm, "start-from")};

        const fs::path nodesFile{input};
        if (!fs::exists(nodesFile)) {
            throw Model::Refused(Model::refuse(RefusalKind::InputMissing, "Nodes input file does not exist."));
        }
        const auto store = given.resume ? given.resume : given.outputDir;
        if (store and !fs::exists(*store)) {
            throw Model::Refused(Model::refuse(RefusalKind::OutputDirMissing, "Output directory does not exist."));
        }
        const auto symptomatic_idp  = loadVectorFromFile<Probability>(symptomatic_idp_path);
        const auto asymptomatic_idp = loadVectorFromFile<Probability>(asymptomatic_idp_path);
        given.inputSha256           = {{Model::input_names::INPUT, sha256HexOfFile(nodesFile)},
                                       {Model::input_names::SYMPTOMATIC_IDP, sha256HexOfFile(symptomatic_idp_path)},
                                       {Model::input_names::ASYMPTOMATIC_IDP, sha256HexOfFile(asymptomatic_idp_path)}};

        // The Run the store holds, and the checkpoint the option given names. Which of them a Run may use is the Run
        // recipe's to rule. A store with no Run has no checkpoint to read: the recipe refuses to resume it.
        const auto stored = store ? Model::readStoredRun(*store) : std::nullopt;
        std::optional<checkpoint::Checkpoint> saved;
        if (given.startFrom) {
            saved = checkpoint::load(*given.startFrom);
        } else if (given.resume and stored) {
            saved = checkpoint::load(*given.resume / checkpoint::LATEST_FILE);
        }

        auto resolved = Model::resolve(given, stored, saved);
        if (auto* refusal = std::get_if<Model::Refusal>(&resolved)) {
            throw Model::Refused(std::move(*refusal));
        }
        auto& recipe        = std::get<Model::Recipe>(resolved);
        const bool resuming = recipe.mode == Model::Mode::Resume;
        Model::ConsoleProgress progress(std::cout);
        for (const auto& notice : recipe.notices) {
            progress.notice(notice);
        }

        if (recipe.manifest.nullModel) {
            progress.notice("Running the null model (no genetics)");
        }

        Model::Inputs inputs{.symptomaticIdp = symptomatic_idp, .asymptomaticIdp = asymptomatic_idp};
        if (nodesFile.extension() == ".gz") {
            std::vector<char> decompressed = decompressGzipFile(nodesFile);
            inputs.input                   = loadJSON(decompressed);
        } else {
            std::ifstream inputFile{nodesFile};
            inputs.input = loadJSON<std::ifstream&>(inputFile);
        }

        const Model::Provenance provenance{.runId = Model::newRunId(), .cliArguments = std::vector<std::string>(argv, argv + argc), .gitCommit = TRANSMISSION_NETWORKS_GIT_COMMIT};
        driving = run.emplace(Model::assemble(recipe, inputs, saved, std::make_unique<ParquetStoreSink>(recipe.store), progress, provenance)).driver.get();
        if (!resuming) {
            progress.notice(fmt::format("Seed Used: {}", run->manifest.seed));
        }

        const auto outcome = run->drive([&]() -> std::optional<std::string> {
            if (!watcher.interrupted()) {
                return std::nullopt;
            }
            return watcher.firstSignal();
        });
        if (std::holds_alternative<Model::Stopped>(outcome)) {
            // The Run stopped cleanly for the first interrupt; report it as a shell reports a process a signal ended.
            return 128 + watcher.firstSignalNumber();
        }

    } catch (const Model::Refused& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return ERROR_IN_COMMAND_LINE;
    } catch (const std::filesystem::filesystem_error& e) {
        failOutput("use", e.path1(), e.code().message());
    } catch (const std::exception& e) {
        std::cerr << "ERROR: unhandled exception: " << e.what() << std::endl;
        return ERROR_UNHANDLED_EXCEPTION;
    }

    return SUCCESS;

}// main
