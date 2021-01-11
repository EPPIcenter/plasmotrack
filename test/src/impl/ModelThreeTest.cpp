//
// Created by Maxwell Murphy on 6/19/20.
//

//#include <boost/filesystem.hpp>
#include <filesystem>
#include <fstream>

#include "gtest/gtest.h"

#include "core/io/parse_json.h"
#include "core/io/path_parsing.h"
#include "core/io/loggers/ValueLogger.h"
#include "core/io/loggers/LambdaLogger.h"

#include "core/samplers/RandomizedScheduler.h"
#include "core/samplers/general/ConstrainedContinuousRandomWalk.h"
#include "core/samplers/general/SALTSampler.h"

#include "core/samplers/genetics/RandomAllelesBitSetSampler.h"

#include "core/samplers/topology/RandomAddEdgeSampler.h"
#include "core/samplers/topology/RandomRemoveEdgeSampler.h"
#include "core/samplers/topology/RandomReverseEdgeSampler.h"
#include "core/samplers/topology/RandomSwapEdgeSampler.h"

#include "impl/model/ModelThree.h"
#include "impl/state/ModelThreeState.h"

namespace fs = std::filesystem;

using namespace transmission_nets::impl;
using namespace transmission_nets::core::io;
using namespace transmission_nets::core::samplers;

TEST(ModelThreeTest, CoreTest) {

    using InfectionEvent = ModelThreeState::InfectionEvent;
    using Locus = ModelThreeState::LocusImpl;
    using GeneticsSampler = genetics::RandomAllelesBitSetSampler<ModelThree, boost::random::mt19937, ModelThreeState::GeneticsImpl>;
    using AddEdgeSampler = topology::RandomAddEdgeSampler<ModelThreeState::MAX_PARENT_SET, ModelThree, boost::random::mt19937, ModelThreeState::InfectionEvent>;
    using RemoveEdgeSampler = topology::RandomRemoveEdgeSampler<ModelThree, boost::random::mt19937, ModelThreeState::InfectionEvent>;
    using ReverseEdgeSampler = topology::RandomReverseEdgeSampler<ModelThreeState::MAX_PARENT_SET, ModelThree, boost::random::mt19937, ModelThreeState::InfectionEvent>;
    using SwapEdgeSampler = topology::RandomSwapEdgeSampler<ModelThree, boost::random::mt19937, ModelThreeState::InfectionEvent>;


    ModelThreeState state;

    const auto testsDir = getPathFromEnvVar("TRANSMISSION_NETWORK_TESTS_DIR");
    const fs::path nodesFile{"resources/JSON/nodes3.json"};
    const fs::path outputDir{"outputs/ModelThreeTests/CoreTest"};

    const fs::path testFilePath = testsDir / nodesFile;
    const fs::path paramOutput = testsDir / outputDir;

    if(!fs::exists(paramOutput)) {
        fs::create_directories(paramOutput);
    }


    if(!fs::exists(testFilePath)) {
        std::cerr << "Nodes test file does not exist." << std::endl;
        exit(1);
    }

    std::ifstream testFile{testFilePath};

    if(!testFile) {
        std::cout << "Cannot open file." << std::endl;
        exit(1);
    }

    auto j = loadJSON(testFile);
    auto loci = parseLociFromJSON<Locus>(j);
    auto infections = parseInfectionsFromJSON<InfectionEvent, Locus>(j, loci);


    state.loci = loci;
    state.infections = infections;

    for(const auto& [locus_label, locus] : state.loci) {
        state.alleleFrequencies.addLocus(locus);
    }

    state.transmissionNetwork.addNodes(state.infections);

    state.observationFalsePositiveRate.initializeValue(.005);
    state.observationFalseNegativeRate.initializeValue(.1);
    state.geometricGenerationProb.initializeValue(.9);
    state.lossProb.initializeValue(.3);
    state.mutationProb.initializeValue(.01);
    state.meanCOI.initializeValue(3);

    ModelThree model(state);

    std::vector<AbstractLogger*> loggers{};
    loggers.push_back(new ValueLogger(paramOutput / "fpr.csv", state.observationFalsePositiveRate));
    loggers.push_back(new ValueLogger(paramOutput / "fnr.csv", state.observationFalseNegativeRate));
    loggers.push_back(new ValueLogger(paramOutput / "geo_gen_prob.csv", state.geometricGenerationProb));
    loggers.push_back(new ValueLogger(paramOutput / "loss_prob.csv", state.lossProb));
    loggers.push_back(new ValueLogger(paramOutput / "mutation_prob.csv", state.mutationProb));
    loggers.push_back(new ValueLogger(paramOutput / "mean_coi.csv", state.meanCOI));
    loggers.push_back(new ValueLogger(paramOutput / "likelihood.csv", model));
    loggers.push_back(new LambdaLogger(paramOutput / "network.csv", [&](){ return state.transmissionNetwork.serialize(); }));
    for(const auto& [locus_label, locus] : state.loci) {
        loggers.push_back(new ValueLogger(paramOutput / (locus_label + "_frequencies.csv"), state.alleleFrequencies.alleleFrequencies(locus)));
    }

    for(const auto& infection : state.infections) {
        for(const auto& [locus_label, locus] : state.loci) {
            if (std::find(infection->loci().begin(), infection->loci().end(), locus) != infection->loci().end()) {
                loggers.push_back(new ValueLogger(paramOutput / "nodes" / (infection->id() + "_" + locus_label + ".csv"), infection->latentGenotype(locus)));
            }
        }
    }

    for(const auto& logger : loggers) {
        logger->clearFile();
    }


    boost::random::mt19937 r;
    RandomizedScheduler scheduler(&r, 500);

    for (int l = 0; l < 10; ++l) {
        scheduler.registerSampler(new AddEdgeSampler(state.transmissionNetwork, model, &r));
        scheduler.registerSampler(new RemoveEdgeSampler(state.transmissionNetwork, model, &r));
        scheduler.registerSampler(new ReverseEdgeSampler(state.transmissionNetwork, model, &r));
        scheduler.registerSampler(new SwapEdgeSampler(state.transmissionNetwork, model, &r));
    }

    for (int k = 0; k < 50000; ++k) {
        scheduler.step();
        if(k % 1000 == 0) {
            std::cout << "(Network) Edge and Genotypes Current LLik: " << model.value() << "\n";
            std::cout << state.transmissionNetwork.serialize() << std::endl;
        }
    }

    for(const auto& [locus_label, locus] : state.loci) {
        scheduler.registerSampler(new SALTSampler<ModelThree>(state.alleleFrequencies.alleleFrequencies(locus), model, &r));
    }

    for(auto &infection : state.infections) {
        for(const auto& [locus_label, locus] : state.loci) {
            if (infection->latentGenotype().contains(locus)) {
                auto &latentGenotype = infection->latentGenotype(locus);
                scheduler.registerSampler(new GeneticsSampler(latentGenotype, model, &r));
            }
        }
    }

//    scheduler.registerSampler(new ConstrainedContinuousRandomWalk(state.observationFalsePositiveRate, model, 0.0, 1.0, &r));
//    scheduler.registerSampler(new ConstrainedContinuousRandomWalk(state.observationFalseNegativeRate, model, 0.0, 1.0, &r));
    scheduler.registerSampler(new ConstrainedContinuousRandomWalk(state.geometricGenerationProb, model, 0.0, 1.0, &r));
    scheduler.registerSampler(new ConstrainedContinuousRandomWalk(state.lossProb, model, 0.0, 1.0, &r));
    scheduler.registerSampler(new ConstrainedContinuousRandomWalk (state.mutationProb, model, 0.0, 0.5, &r));
    scheduler.registerSampler(new ConstrainedContinuousRandomWalk(state.meanCOI, model, 0.0, std::numeric_limits<double>::infinity(), &r, .1, .01, 1));




    for (int k = 0; k < 5000; ++k) {
        scheduler.step();
        if(k % 100 == 0) {
            std::cout << "Edge and Genotypes Current LLik: " << model.value() << "\n";
            std::cout << state.transmissionNetwork.serialize() << std::endl;
            std::cout << state.observationFalsePositiveRate.value() << " ";
            std::cout << state.observationFalseNegativeRate.value() << " ";
            std::cout << state.lossProb.value() << " ";
            std::cout << state.mutationProb.value() << std::endl;
        }
    }

    for (int i = 0; i < 5000; ++i) {
        scheduler.step();
        if (i % 100 == 0) {
            std::cout << "(Writing) Current LLik: " << model.value() <<  std::endl;
            for (const auto& logger : loggers) {
                logger->logValue();
            }
        }
    }

    std::cout << "Completed" << "\n";
}