//
// Created by Maxwell Murphy on 6/5/20.
//

//#include <boost/filesystem.hpp>

#include <filesystem>
#include <random>

#include "gtest/gtest.h"

#include "core/io/parse_json.h"
#include "core/io/loggers/ValueLogger.h"
#include "core/io/loggers/ParentSetDistLogger.h"
#include "core/io/path_parsing.h"
#include "core/io/serialize.h"


#include "core/samplers/RandomizedScheduler.h"
#include "core/samplers/general/ConstrainedContinuousRandomWalk.h"
#include "core/samplers/general/SALTSampler.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler.h"
#include "core/samplers/genetics/ZanellaAllelesBitSetSampler.h"
#include "core/samplers/topology/OrderSampler.h"
#include "core/samplers/topology/ZanellaNeighborOrderSampler.h"
#include "core/samplers/topology/ZanellaOrderSampler.h"

#include "impl/model/ModelTwo.h"


namespace fs = std::filesystem;

using namespace transmission_nets::impl;
using namespace transmission_nets::core::io;
using namespace transmission_nets::core::samplers;

TEST(ModelTwoTest, CoreTest) {

    using InfectionEvent = ModelTwo::InfectionEvent;
    using Locus = ModelTwo::LocusImpl;
//    using InformedGeneticsSampler = genetics::ZanellaAllelesBitSetSampler<ModelTwo, boost::random::mt19937, ModelTwo::GeneticsImpl>;
//    using GeneticsSampler = genetics::RandomAllelesBitSetSampler<ModelTwo, boost::random::mt19937, ModelTwo::GeneticsImpl>;

    const auto testsDir = getPathFromEnvVar("TRANSMISSION_NETWORK_TESTS_DIR");
    const fs::path nodesFile{"resources/JSON/nodes4.json"};
    const fs::path outputDir{"outputs/ModelTwoTests/CoreTest"};

    const fs::path testFilePath = testsDir / nodesFile;
    const fs::path paramOutput = testsDir / outputDir / "params";
    const fs::path statOutput = testsDir / outputDir / "stats";

    if(!fs::exists(paramOutput)) {
        fs::create_directories(paramOutput);
    }

    if(!fs::exists(statOutput)) {
        fs::create_directories(statOutput);
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
    auto disallowedParents = parseDisallowedParentsFromJSON(j, infections);


    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    auto rng = std::default_random_engine {seed};
    std::shuffle(infections.begin(), infections.end(), rng);

    ModelTwo model(loci, infections, disallowedParents);

    std::vector<AbstractLogger*> loggers{};
    loggers.push_back(new ValueLogger(paramOutput / "fpr.csv", model.state.observationFalsePositiveRate));
    loggers.push_back(new ValueLogger(paramOutput / "fnr.csv", model.state.observationFalseNegativeRate));
    loggers.push_back(new ValueLogger(paramOutput / "geo_gen_prob.csv", model.state.geometricGenerationProb));
    loggers.push_back(new ValueLogger(paramOutput / "loss_prob.csv", model.state.lossProb));
    loggers.push_back(new ValueLogger(paramOutput / "mutation_prob.csv", model.state.mutationProb));
    loggers.push_back(new ValueLogger(paramOutput / "mean_coi.csv", model.state.meanCOI));
    loggers.push_back(new ValueLogger(paramOutput / "infection_order.csv", model.state.infectionEventOrdering));
    loggers.push_back(new ValueLogger(paramOutput / "likelihood.csv", model));
    for(const auto& [locus_label, locus] : model.state.loci) {
        loggers.push_back(new ValueLogger(paramOutput / (locus_label + "_frequencies.csv"), model.state.alleleFrequencies.alleleFrequencies(locus)));
    }

    for(const auto& infection : model.state.infections) {
        for(const auto& [locus_label, locus] : model.state.loci) {
            if (std::find(infection->loci().begin(), infection->loci().end(), locus) != infection->loci().end()) {
                loggers.push_back(new ValueLogger(paramOutput / "nodes" / (infection->id() + "_" + locus_label + ".csv"), infection->latentGenotype(locus)));
            }
        }
    }

    for (const auto& tp : model.transmissionProcessList) {
        loggers.push_back(new ParentSetDistLogger(statOutput / (tp->child_.id() + "_ps.csv"), tp));
    }

    for(const auto& logger : loggers) {
        logger->clearFile();
    }

    boost::random::mt19937 r;
    RandomizedScheduler scheduler(&r, 5000);
//
    scheduler.registerSampler({.sampler = new ConstrainedContinuousRandomWalk(model.state.observationFalsePositiveRate, model, 0.0, 0.1, &r, .01),
//                               .adaptationStart = 1000,
//                               .adaptationEnd = 3000,
                               .updateStart = 1000,
                               .weight = 20});

    scheduler.registerSampler({.sampler = new ConstrainedContinuousRandomWalk(model.state.observationFalseNegativeRate, model, 0.0, 0.2, &r, .01),
//                               .adaptationStart = 1000,
//                               .adaptationEnd = 3000,
                               .updateStart = 1000,
                               .weight = 20});
//
    scheduler.registerSampler({.sampler = new ConstrainedContinuousRandomWalk(model.state.geometricGenerationProb, model, 0.0, 1.0, &r, .01, .01 ,2),
//                               .adaptationStart = 1000,
//                               .adaptationEnd = 3000,
//                               .scaledAdaptation = true,
                               .weight = 20});

    scheduler.registerSampler({.sampler = new ConstrainedContinuousRandomWalk(model.state.lossProb, model, 0.0, 1.0, &r, .01, .01, 2),
//                               .adaptationStart = 1000,
//                               .adaptationEnd = 3000,
//                               .updateStart = 100,
//                               .scaledAdaptation = true,
                               .weight = 20});

//    scheduler.registerSampler({.sampler = new ProbabilitySampler(model.state.mutationProb, model, 0.0, .1, &r, .1, .01, 1),
//                               .adaptationStart = 1000,
//                               .adaptationEnd = 3000,
//                               .weight = 5});

    scheduler.registerSampler({.sampler = new ConstrainedContinuousRandomWalk(model.state.meanCOI, model, 0.0, double(model.MAX_COI), &r, .1, .01, 1),
//                               .adaptationStart = 1000,
//                               .adaptationEnd = 3000,
//                               .scaledAdaptation = true,
                               .weight = 20});

//    scheduler.registerSampler({.sampler = new ZanellaOrderSampler(model.state.infectionEventOrdering, model, &r, 5),
////                               .updateStart = 0,
//                               .weight = 50});
////
//    scheduler.registerSampler({.sampler = new ZanellaNeighborOrderSampler(model.state.infectionEventOrdering, model, &r),
////                               .updateStart = 0,
//                               .weight = 1});

    scheduler.registerSampler({.sampler = new topology::OrderSampler(model.state.infectionEventOrdering, model, &r, 50),
//                           .adaptationStart = 1000,
//                           .adaptationEnd = 3000,
                           .weight = 500});

//    for (int l = 1; l < (int) model.state.infections.size() / 2; ++l) {
//        scheduler.registerSampler({.sampler = new OrderSampler(model.state.infectionEventOrdering, model, &r, 5),
//                                   .adaptationStart = 1000,
//                                   .adaptationEnd = 3000,
//                                   .weight = 400});
//    }

    for (auto &infection : model.state.infections) {
        for (const auto &[locus_label, locus] : model.state.loci) {
            if (infection->latentGenotype().contains(locus)) {
                auto &latentGenotype = infection->latentGenotype(locus);
                scheduler.registerSampler({.sampler = new genetics::ZanellaAllelesBitSetSampler(latentGenotype, model, &r),
//                                           .updateStart = 0,
//                                           .updateEnd = 1000,
                                           .weight = .02});
                scheduler.registerSampler({.sampler = new genetics::RandomAllelesBitSetSampler(latentGenotype, model, &r),
//                                           .updateStart = 0,
                                           .weight = 2});
            }
        }
    }

    for (const auto &[locus_label, locus] : model.state.loci) {
        scheduler.registerSampler({.sampler = new SALTSampler<ModelTwo>(model.state.alleleFrequencies.alleleFrequencies(locus), model, &r, .05, .01, .5),
//                                   .adaptationStart = 100,
//                                   .adaptationEnd = 3000,
                                   .weight = 1});
    }

    std::cout << "Current LLik: " << model.value() << std::endl;
    for (int k = 0; k < 30000; ++k) {
        scheduler.step();
        std::cout << "(k=" << k <<  ") Current LLik: " << model.value() << std::endl;
        if(k % 10 == 0) {
            std::cout << "Current LLik: " << model.value() << std::endl;
            std::cout << serialize(model.state.infectionEventOrdering.value()) << std::endl;
            std::cout << serialize(model.state.alleleFrequencies.alleleFrequencies(model.state.loci.begin()->second).value()) << std::endl;
            std::cout << model.state.meanCOI.value() << " ";
            std::cout << model.state.geometricGenerationProb.value() << " ";
            std::cout << model.state.observationFalsePositiveRate.value() << " ";
            std::cout << model.state.observationFalseNegativeRate.value() << " ";
            std::cout << model.state.lossProb.value() << std::endl;
            std::cout << model.state.mutationProb.value() << std::endl;
            std::cout << model.transmissionProcessList[0]->child_.id() << ": ";
            auto distr = model.transmissionProcessList[0]->calcParentSetDist();
            for (const auto& [llik, ps] : distr.parentSetLliks) {
                std::cout << serialize(ps) << ": " << std::exp(llik - distr.totalLlik) << std::endl;
            }
            std::cout << "{S}: " << std::exp(distr.sourceLlik - distr.totalLlik) << std::endl;
        }
    }

    for (int i = 0; i < 50000; ++i) {
        scheduler.step();
        std::cout << "(i=" << i <<  ") Current LLik: " << model.value() << std::endl;
        if (i % 10 == 0) {
            std::cout << "Current LLik: " << model.value() << std::endl;
            std::cout << serialize(model.state.infectionEventOrdering.value()) << std::endl;
            std::cout << serialize(model.state.alleleFrequencies.alleleFrequencies(model.state.loci.begin()->second).value()) << std::endl;
            std::cout << model.state.meanCOI.value() << " ";
            std::cout << model.state.geometricGenerationProb.value() << " ";
            std::cout << model.state.observationFalsePositiveRate.value() << " ";
            std::cout << model.state.observationFalseNegativeRate.value() << " ";
            std::cout << model.state.lossProb.value() << " ";
            std::cout << model.state.mutationProb.value() << std::endl;
            std::cout << model.transmissionProcessList[0]->child_.id() << ": ";
            auto distr = model.transmissionProcessList[0]->calcParentSetDist();
            for (const auto& [llik, ps] : distr.parentSetLliks) {
                std::cout << serialize(ps) << ": " << std::exp(llik - distr.totalLlik) << std::endl;
            }
            std::cout << "{S}: " << std::exp(distr.sourceLlik - distr.totalLlik) << std::endl;

            for (const auto& logger : loggers) {
                logger->logValue();
            }
        }
    }

    std::cout << "Completed" << std::endl;
}