//
// Created by mmurphy on 11/1/21.
//

#include "core/io/parse_json.h"
#include "core/samplers/meta/ReplicaExchange.h"
#include "core/utils/timers.h"
#include "impl/model/ModelEight/Model.h"
#include "impl/model/ModelEight/ModelLogger.h"
#include "impl/model/ModelEight/SampleScheduler.h"
#include "impl/model/ModelEight/State.h"
#include "impl/model/ModelEight/StateLogger.h"

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
using namespace transmission_nets::core::utils;

namespace {
    const size_t SUCCESS                   = 0;
    const size_t ERROR_IN_COMMAND_LINE     = 1;
    const size_t ERROR_UNHANDLED_EXCEPTION = 2;

}// namespace

// Global interrupt tracker
bool interrupted = false;

// Global replica exchange object
std::unique_ptr<ReplicaExchange<ModelEight::State, ModelEight::Model, ModelEight::SampleScheduler, ModelEight::ModelLogger, ModelEight::StateLogger>> repex;

void finalize_output(int signal_num) {
    if (interrupted) {
        fmt::print("Killing process, output inconsistent...\n");
        exit(SUCCESS);
    }
    interrupted = true;
    fmt::print("Interrupt signal {} received. Finalizing output...", signal_num);
}

//using Time = std::chrono::high_resolution_clock;
//using dsec = std::chrono::duration<double>;

int main(int argc, char** argv) {
    signal(SIGINT, finalize_output);
    signal(SIGQUIT, finalize_output);
    signal(SIGABRT, finalize_output);

//    try {
        int num_chains;
        unsigned int num_cores;
        double gradient;

        int burnin;
        int sample;
        int thin;
        long seed;
        bool null_model;
        std::string input;
        std::string output_dir;

        namespace po = boost::program_options;
        po::options_description desc("Options");
        auto opts = desc.add_options();
        opts("help", "Runs the model six implementation");
        opts("burnin,b", po::value<int>(&burnin)->default_value(5000), "Number of steps to be used for burnin");
        opts("sample,s", po::value<int>(&sample)->default_value(10000), "Total number of steps to be used for sampling");
        opts("thin,t", po::value<int>(&thin)->default_value(1000), "Number of steps to be thinned");
        opts("numchains,n", po::value<int>(&num_chains)->default_value(1), "Number of chains to run in replica exchange algorithm.");
        opts("numcores,c", po::value<unsigned int>(&num_cores)->default_value(1), "Number of cores to use in replica exchange algorithm.");
        opts("gradient,g", po::value<double>(&gradient)->default_value(1), "Temperature gradient to use in replica exchange algorithm");
        opts("seed", po::value<long>(&seed)->default_value(-1), "Seed used in random number generator. Note that if numchains > 1 then there is no guarantee of reproducibility. A value of -1 indicates generate a random seed.");
        opts("hotload,h", "Hotload parameters from the output directory");
        opts("input,i", po::value<std::string>(&input)->required(), "Input file");
        opts("output-dir,o", po::value<std::string>(&output_dir)->required(), "Output directory");
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
                std::cout << "Model Eight Implementation" << std::endl
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
            std::cerr << "Nodes input file does not exist." << std::endl;
            return ERROR_IN_COMMAND_LINE;
        }

        if (!fs::exists(outputDir)) {
            std::cerr << "Output directory does not exist." << std::endl;
            return ERROR_IN_COMMAND_LINE;
        }

        if (null_model) {
            fmt::print("Running the null model (no genetics)\n");
        }

        std::ifstream inputFile{nodesFile};
        auto j = loadJSON(inputFile);

        if (seed == -1) {
            seed = timers::time().time_since_epoch().count();
        }

        fmt::print("Seed Used: {}\n", seed);
        auto r = std::make_shared<boost::random::mt19937>(seed);
        repex  = std::make_unique<ReplicaExchange<ModelEight::State, ModelEight::Model, ModelEight::SampleScheduler, ModelEight::ModelLogger, ModelEight::StateLogger>>(num_chains, thin, gradient, r, outputDir, hotload, null_model, num_cores, j);

        repex->logState();
        repex->finalize();

        double totalSamples = 0;
        double totalSeconds = 0;
        double samplesPerSecond = 0;

        fmt::print("Starting Llik: {0:.2f}\n", repex->hotValue());
        for (int kk = 0; kk < burnin; ++kk) {
            if (interrupted) {
                repex->logModel();
                repex->logState();
                break;
            }

            auto t0 = timers::time();
            repex->sample();
            auto t1 = timers::time();

            timers::dsec ds = t1 - t0;
            totalSamples += thin;
            totalSeconds += ds.count();
            samplesPerSecond = (totalSamples / totalSeconds);

            fmt::print("(b={0})", kk);
            repex->printModelLlik();
            fmt::print("({0:.2f} samples/sec)\n", samplesPerSecond);

//            fmt::print("(b={0}) Current Llik: {1:.2f} ({2} samples/sec)\n", kk, repex->hotValue(), samplesPerSecond);
        }


        for (int jj = 0; jj < sample; ++jj) {
            if (interrupted) {
                break;
            }

            auto t0 = timers::time();
            repex->sample();
            auto t1 = timers::time();

            timers::dsec ds = t1 - t0;
            totalSamples += thin;
            totalSeconds += ds.count();
            samplesPerSecond = thin / ds.count();

            repex->logModel();
            repex->logState();

            fmt::print("(s={0}) ", jj);
            repex->printModelLlik();
            fmt::print(" ({0:.2f} samples/sec)\n", samplesPerSecond);

//            fmt::print("(s={0}) Current Llik: {1:.2f} ({2} samples/sec)\n", jj, repex->hotValue(), samplesPerSecond);
        }

        repex->finalize();

//    } catch (std::exception& e) {
//        std::cerr << "Unhandled Exception reached the top of main: "
//                  << e.what() << ", application will now exit" << std::endl;
//        return ERROR_UNHANDLED_EXCEPTION;
//    }

    return SUCCESS;

}// main