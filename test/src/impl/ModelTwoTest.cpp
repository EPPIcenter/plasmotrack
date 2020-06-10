//
// Created by Maxwell Murphy on 6/5/20.
//

#include <boost/filesystem.hpp>

#include "gtest/gtest.h"

#include "core/utils/io/parse_json.h"
#include "core/utils/io/Loggers/ValueLogger.h"
#include "core/utils/io/path_parsing.h"

#include "core/samplers/ConstrainedContinuousRandomWalk.h"
#include "core/samplers/SALTSampler.h"
#include "core/samplers/genetics/RandomAllelesBitSetSampler.h"
#include "core/samplers/OrderSampler.h"
#include "core/samplers/RandomizedScheduler.h"

#include "impl/model/ModelTwo.h"
#include "impl/state/ModelTwoState.h"


namespace fs = boost::filesystem;

TEST(ModelTwoTest, CoreTest) {

    using InfectionEvent = ModelTwoState::InfectionEvent;
    using Locus = ModelTwoState::LocusImpl;
    using ZeroOneSampler = ConstrainedContinuousRandomWalk<0, 1, ModelTwo, boost::random::mt19937>;
    using ZeroBoundedSampler = ConstrainedContinuousRandomWalk<0, std::numeric_limits<int>::max(), ModelTwo, boost::random::mt19937>;
    using GeneticsSampler = RandomAllelesBitSetSampler<ModelTwo, boost::random::mt19937, ModelTwoState::GeneticsImpl>;

    ModelTwoState state;

    const auto testsDir = getPathFromEnvVar("TRANSMISSION_NETWORK_TESTS_DIR");
    const fs::path nodesFile{"resources/JSON/nodes.json"};
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

    fs::ifstream testFile{testFilePath};

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

    state.observationFalsePositiveRate.initializeValue(.15);
    state.observationFalseNegativeRate.initializeValue(.15);
    state.geometricGenerationProb.initializeValue(.5);
    state.lossProb.initializeValue(.001);
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
    RandomizedScheduler scheduler(&r);
    for (int l = 1; l < (int)state.infections.size() / 2; ++l) {
        scheduler.registerSampler(new OrderSampler(state.infectionEventOrdering, model, &r, l));
    }
    scheduler.registerSampler(new OrderSampler(state.infectionEventOrdering, model, &r, 1));
    scheduler.registerSampler(new OrderSampler(state.infectionEventOrdering, model, &r, 2));
    scheduler.registerSampler(new OrderSampler(state.infectionEventOrdering, model, &r, 3));
    scheduler.registerSampler(new OrderSampler(state.infectionEventOrdering, model, &r, 4));
    scheduler.registerSampler(new OrderSampler(state.infectionEventOrdering, model, &r, 5));
    scheduler.registerSampler(new OrderSampler(state.infectionEventOrdering, model, &r, 6));
    scheduler.registerSampler(new OrderSampler(state.infectionEventOrdering, model, &r, 7));
    scheduler.registerSampler(new OrderSampler(state.infectionEventOrdering, model, &r, 8));
    scheduler.registerSampler(new OrderSampler(state.infectionEventOrdering, model, &r, 9));
//    scheduler.registerSampler(new ZeroOneSampler(state.observationFalsePositiveRate, model, &r));
//    scheduler.registerSampler(new ZeroOneSampler(state.observationFalseNegativeRate, model, &r));
//    scheduler.registerSampler(new ZeroOneSampler(state.geometricGenerationProb, model, &r));
//    scheduler.registerSampler(new ZeroOneSampler(state.lossProb, model, &r));
//    scheduler.registerSampler(new ZeroOneSampler (state.mutationProb, model, &r));
//    scheduler.registerSampler(new ZeroBoundedSampler(state.meanCOI, model, &r));


    for(auto &infection : state.infections) {
        for(const auto& [locus_label, locus] : state.loci) {
            if (infection->latentGenotype().contains(locus)) {
                auto &latentGenotype = infection->latentGenotype(locus);
                scheduler.registerSampler(new GeneticsSampler(latentGenotype, model, &r));
            }
        }
    }

//    const auto genotypeSampler = dynamic_cast<GeneticsSampler*>(scheduler.samplers().back());

    for(const auto& [locus_label, locus] : state.loci) {
        scheduler.registerSampler(new SALTSampler<ModelTwo>(state.alleleFrequencies.alleleFrequencies(locus), model, &r));
    }

    for (int k = 0; k < 50000; ++k) {
//        scheduler.update();
        scheduler.updateAndAdapt();
        if(k % 10 == 0) {
            std::cout << "Current LLik: " << model.value() << std::endl;
//            std::cout << "Genotype Sampler: " << genotypeSampler->acceptanceRate() << std::endl;
        }
    }

    for (int i = 0; i < 50000; ++i) {
        scheduler.update();
        if (i % 10 == 0) {
            for (const auto& logger : loggers) {
                logger->logValue();
            }
        }
    }

    std::cout << "Completed" << std::endl;
}