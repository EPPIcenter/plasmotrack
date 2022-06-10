//
// Created by mmurphy on 11/1/21.
//

#include "core/io/parse_json.h"
#include "core/samplers/meta/ReplicaExchange.h"
#include "impl/model/ModelSeven/Model.h"
#include "impl/model/ModelSeven/ModelLogger.h"
#include "impl/model/ModelSeven/SampleScheduler.h"
#include "impl/model/ModelSeven/State.h"
#include "impl/model/ModelSeven/StateLogger.h"

#include <boost/program_options.hpp>
#include <fmt/core.h>

#include <csignal>
#include <filesystem>
#include <iostream>
#include <string>


using namespace transmission_nets::impl;
using namespace transmission_nets::core::io;
using namespace transmission_nets::core::computation;
using namespace transmission_nets::core::samplers;

namespace {
    const size_t SUCCESS                   = 0;
    const size_t ERROR_IN_COMMAND_LINE     = 1;
    const size_t ERROR_UNHANDLED_EXCEPTION = 2;

}// namespace

// Global interrupt tracker
bool interrupted = false;

// Global replica exchange object
std::unique_ptr<ReplicaExchange<ModelSeven::State, ModelSeven::Model, ModelSeven::SampleScheduler, ModelSeven::ModelLogger, ModelSeven::StateLogger>> repex;

void finalize_output(int signal_num) {
    interrupted = true;
    fmt::print("Interrupt signal {} received. Finalizing output...", signal_num);
}

int main(int argc, char** argv) {
    signal(SIGINT, finalize_output);
    signal(SIGQUIT, finalize_output);
    signal(SIGABRT, finalize_output);

    try {
        int numChains;
        double gradient;

        int burnin;
        int sample;
        int thin;
        int seed;
        std::string input;
        std::string output_dir;

        namespace po = boost::program_options;
        po::options_description desc("Options");
        auto opts = desc.add_options();
        opts("help", "Runs the model six implementation");
        opts("burnin,b", po::value<int>(&burnin)->default_value(5000), "Number of steps to be used for burnin");
        opts("sample,s", po::value<int>(&sample)->default_value(10000), "Total number of steps to be used for sampling");
        opts("thin,t", po::value<int>(&thin)->default_value(1000), "Number of steps to be thinned");
        opts("numchains,n", po::value<int>(&numChains)->default_value(1), "Number of chains to run in replica exchange algorithm. Do not exceed the number of threads available.");
        opts("gradient,g", po::value<double>(&gradient)->default_value(1), "Temperature gradient to use in replica exchange algorithm");
        opts("seed", po::value<int>(&seed)->default_value(-1), "Seed used in random number generator. Note that if numchains > 1 then there is no guarantee of reproducibility. A value of -1 indicates generate a random seed.");
        opts("hotload,h", "Hotload parameters from the output directory");
        opts("input,i", po::value<std::string>(&input)->required(), "Input file");
        opts("output-dir,o", po::value<std::string>(&output_dir)->required(), "Output directory");


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
                std::cout << "Model Seven Implementation" << std::endl
                          << desc << std::endl;
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

        bool hotload = false;
        if (vm.count("hotload")) {
            hotload = true;
        }

        namespace fs = std::filesystem;

        const fs::path nodesFile{input};
        const fs::path outputDir{output_dir};

        if (!fs::exists(nodesFile)) {
            std::cerr << "Nodes test file does not exist." << std::endl;
            return ERROR_IN_COMMAND_LINE;
        }

        if (!fs::exists(outputDir)) {
            std::cerr << "Output directory does not exist." << std::endl;
            return ERROR_IN_COMMAND_LINE;
        }

        std::ifstream inputFile{nodesFile};
        auto j = loadJSON(inputFile);

        if (seed == -1) {
            seed = std::chrono::system_clock::now().time_since_epoch().count();
        }

        fmt::print("Seed Used: {}\n", seed);
        auto r = std::make_shared<boost::random::mt19937>(seed);
        repex  = std::make_unique<ReplicaExchange<ModelSeven::State, ModelSeven::Model, ModelSeven::SampleScheduler, ModelSeven::ModelLogger, ModelSeven::StateLogger>>(numChains, thin, gradient, r, outputDir, hotload, j);

        repex->logState();
        repex->finalize();
        fmt::print("Starting Llik: {0:.2f}\n", repex->hotValue());
        for (int kk = 0; kk < burnin; ++kk) {
            if (interrupted) {
                repex->logModel();
                repex->logState();
                break;
            }
            repex->sample();
            fmt::print("(b={0}) Current Llik: {1:.2f}\n", kk, repex->hotValue());
        }


        for (int jj = 0; jj < sample; ++jj) {
            if (interrupted) {
                break;
            }
            repex->sample();
            repex->logModel();
            repex->logState();
            fmt::print("(s={0}) Current Llik: {1:.2f}\n", jj, repex->hotValue());
        }

        repex->finalize();

    } catch (std::exception& e) {
        std::cerr << "Unhandled Exception reached the top of main: "
                  << e.what() << ", application will now exit" << std::endl;
        return ERROR_UNHANDLED_EXCEPTION;
    }

    return SUCCESS;

}// main