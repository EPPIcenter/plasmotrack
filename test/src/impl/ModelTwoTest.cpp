//
// Created by Maxwell Murphy on 6/5/20.
//

//#include <boost/filesystem.hpp>

#include <filesystem>
#include <fstream>

#include "gtest/gtest.h"

#include "core/io/parse_json.h"
#include "core/io/loggers/ValueLogger.h"
#include "core/io/path_parsing.h"
#include "core/io/serialize.h"


#include "core/samplers/DoubleConstrainedContinuousRandomWalk.h"
#include "core/samplers/ConstrainedContinuousRandomWalk.h"
#include "core/samplers/SALTSampler.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler.h"
#include "core/samplers/OrderSampler.h"
#include "core/samplers/RandomizedScheduler.h"

#include "impl/model/ModelTwo.h"
#include "impl/state/ModelTwoState.h"



namespace fs = std::filesystem;

using namespace transmission_nets::impl;
using namespace transmission_nets::core::io;
using namespace transmission_nets::core::samplers;

TEST(ModelTwoTest, CoreTest) {

    using InfectionEvent = ModelTwoState::InfectionEvent;
    using Locus = ModelTwoState::LocusImpl;
    using ProbabilitySampler = DoubleConstrainedContinuousRandomWalk<ModelTwo, boost::random::mt19937>;
    using ZeroBoundedSampler = ConstrainedContinuousRandomWalk<0, std::numeric_limits<int>::max(), ModelTwo, boost::random::mt19937>;
    using GeneticsSampler = genetics::RandomAllelesBitSetSampler<ModelTwo, boost::random::mt19937, ModelTwoState::GeneticsImpl>;

    ModelTwoState state;

    const auto testsDir = getPathFromEnvVar("TRANSMISSION_NETWORK_TESTS_DIR");
    const fs::path nodesFile{"resources/JSON/nodes3.json"};
    const fs::path outputDir{"outputs/ModelTwoTests/CoreTest"};

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

    state.infectionEventOrdering.addElements(state.infections);

    state.observationFalsePositiveRate.initializeValue(.001);
    state.observationFalseNegativeRate.initializeValue(.1);
    state.geometricGenerationProb.initializeValue(.8);
    state.lossProb.initializeValue(.3);
    state.mutationProb.initializeValue(.0005);
    state.meanCOI.initializeValue(1);

    ModelTwo model(state);

    std::vector<AbstractLogger*> loggers{};
    loggers.push_back(new ValueLogger(paramOutput / "fpr.csv", state.observationFalsePositiveRate));
    loggers.push_back(new ValueLogger(paramOutput / "fnr.csv", state.observationFalseNegativeRate));
    loggers.push_back(new ValueLogger(paramOutput / "geo_gen_prob.csv", state.geometricGenerationProb));
    loggers.push_back(new ValueLogger(paramOutput / "loss_prob.csv", state.lossProb));
    loggers.push_back(new ValueLogger(paramOutput / "mutation_prob.csv", state.mutationProb));
    loggers.push_back(new ValueLogger(paramOutput / "mean_coi.csv", state.meanCOI));
    loggers.push_back(new ValueLogger(paramOutput / "infection_order.csv", state.infectionEventOrdering));
    loggers.push_back(new ValueLogger(paramOutput / "likelihood.csv", model));
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
    RandomizedScheduler scheduler(&r, 2000);

    scheduler.registerSampler({.sampler = new ProbabilitySampler(state.observationFalsePositiveRate, model, 0.0, .2, &r),
                               .adaptationStart = 0,
                               .adaptationEnd = 0,
                               .weight = 1});

    scheduler.registerSampler({.sampler = new ProbabilitySampler(state.observationFalseNegativeRate, model, 0.0, .2, &r),
                               .adaptationStart = 0,
                               .adaptationEnd = 0,
                               .weight = 1});

    scheduler.registerSampler({.sampler = new ProbabilitySampler(state.geometricGenerationProb, model, 0.0, 1.0, &r),
                               .adaptationStart = 0,
                               .adaptationEnd = 0,
                               .weight = 1});

    scheduler.registerSampler({.sampler = new ProbabilitySampler(state.lossProb, model, 0.0, 1.0, &r),
                               .adaptationStart = 0,
                               .adaptationEnd = 0,
                               .weight = 1});

    scheduler.registerSampler({.sampler = new ProbabilitySampler(state.mutationProb, model, 0.0, .1, &r),
                               .adaptationStart = 0,
                               .adaptationEnd = 0,
                               .weight = 1});

    scheduler.registerSampler({.sampler = new ZeroBoundedSampler(state.meanCOI, model, &r),
                               .adaptationStart = 0,
                               .adaptationEnd = 0,
                               .weight = 1});

    for (int l = 1; l < (int) state.infections.size() / 2; ++l) {
        scheduler.registerSampler({.sampler = new OrderSampler(state.infectionEventOrdering, model, &r, l),
                                   .adaptationStart = 0,
                                   .adaptationEnd = 0,
                                   .weight = 1});
    }

    for (auto &infection : state.infections) {
        for (const auto &[locus_label, locus] : state.loci) {
            if (infection->latentGenotype().contains(locus)) {
                auto &latentGenotype = infection->latentGenotype(locus);
                scheduler.registerSampler({
                        .sampler = new GeneticsSampler(latentGenotype, model, &r),
                        .updateStart = 0,
                        .weight = 1
                });
            }
        }
    }

    for(const auto& [locus_label, locus] : state.loci) {
        scheduler.registerSampler({.sampler = new SALTSampler<ModelTwo>(state.alleleFrequencies.alleleFrequencies(locus), model, &r),
            .adaptationStart = 0,
            .adaptationEnd = 20000,
            .scaledAdaptation = true,
            .weight = 1});
    }

    std::cout << "Current LLik: " << model.value() << std::endl;
    for (int k = 0; k < 20000; ++k) {
        scheduler.step();
        std::cout << "Current LLik: " << model.value() << std::endl;
        if(k % 10 == 0) {
            std::cout << "Current LLik: " << model.value() << std::endl;
            std::cout << serialize(state.infectionEventOrdering.value()) << std::endl;
            std::cout << serialize(state.alleleFrequencies.alleleFrequencies(state.loci.begin()->second).value()) << std::endl;
            std::cout << state.observationFalsePositiveRate.value() << " ";
            std::cout << state.observationFalseNegativeRate.value() << " ";
            std::cout << state.lossProb.value() << " ";
            std::cout << state.mutationProb.value() << std::endl;
        }
    }

    for (int i = 0; i < 50000; ++i) {
        scheduler.step();
        std::cout << "Current LLik: " << model.value() << std::endl;
        if (i % 10 == 0) {
            std::cout << "Current LLik: " << model.value() << std::endl;
            std::cout << serialize(state.infectionEventOrdering.value()) << std::endl;
            std::cout << state.observationFalsePositiveRate.value() << " ";
            std::cout << state.observationFalseNegativeRate.value() << " ";
            std::cout << state.lossProb.value() << " ";
            std::cout << state.mutationProb.value() << std::endl;
            for (const auto& logger : loggers) {
                logger->logValue();
            }
        }
    }

    std::cout << "Completed" << std::endl;
}