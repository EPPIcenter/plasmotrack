////
//// Created by Maxwell Murphy on 4/20/20.
////
//
////#include <boost/filesystem.hpp>
//#include <filesystem>
//#include <fstream>
//
//#include "gtest/gtest.h"
//
//#include "core/io/parse_json.h"
//#include "core/io/loggers/ValueLogger.h"
//#include "core/io/path_parsing.h"
//
//#include "core/samplers/RandomizedScheduler.h"
//#include "core/samplers/general/ConstrainedContinuousRandomWalk.h"
//#include "core/samplers/general/SALTSampler.h"
//#include "core/samplers/genetics/RandomAllelesBitSetSampler.h"
//#include "core/samplers/topology/OrderSampler.h"
//
//#include "impl/model/ModelOne.h"
//#include "impl/state/ModelOneState.h"
//
//using namespace transmission_nets::impl;
//using namespace transmission_nets::core::io;
//using namespace transmission_nets::core::samplers;
//using namespace transmission_nets::core::samplers::topology;
////using namespace transmission_nets::core::samplers::genetics;
//
//namespace fs = std::filesystem;
//
//TEST(ModelOneTest, CoreTest) {
//
//    using InfectionEvent = ModelOneState::InfectionEvent;
//    using Locus = ModelOneState::LocusImpl;
//    using GeneticsSampler = genetics::RandomAllelesBitSetSampler<ModelOne, boost::random::mt19937, ModelOneState::GeneticsImpl>;
//
//    ModelOneState state;
//
//    const auto testsDir = getPathFromEnvVar("TRANSMISSION_NETWORK_TESTS_DIR");
//    const fs::path nodesFile{"resources/JSON/nodes.json"};
//    const fs::path outputDir{"outputs/ModelOneTests/CoreTest"};
//
//    const fs::path testFilePath = testsDir / nodesFile;
//    const fs::path paramOutput = testsDir / outputDir;
//
//    if(!fs::exists(paramOutput)) {
//        fs::create_directories(paramOutput);
//    }
//
//
//    if(!fs::exists(testFilePath)) {
//        std::cerr << "Nodes test file does not exist." << std::endl;
//        exit(1);
//    }
//
//    std::ifstream testFile{testFilePath};
//
//    if(!testFile) {
//        std::cout << "Cannot open file." << std::endl;
//        exit(1);
//    }
//
//    auto j = loadJSON(testFile);
//    auto loci = parseLociFromJSON<Locus>(j);
//    auto infections = parseInfectionsFromJSON<InfectionEvent, Locus>(j, loci);
//
//
//    state.loci = loci;
//    state.infections = infections;
//
//    for(const auto& [locus_label, locus] : state.loci) {
//        state.alleleFrequencies.addLocus(locus);
//    }
//
//    state.infectionEventOrdering.addElements(state.infections);
//
//    state.observationFalsePositiveRate.initializeValue(.025);
//    state.observationFalseNegativeRate.initializeValue(.025);
//    state.geometricGenerationProb.initializeValue(.5);
//    state.ztMultiplicativeBinomialProb.initializeValue(.9);
//    state.ztMultiplicativeBinomialAssoc.initializeValue(1.0);
//    state.geometricCOIProb.initializeValue(.9);
//
//    ModelOne model(state);
//
//    std::vector<FileOutput *> loggers{};
//    loggers.push_back(new ValueLogger(paramOutput / "fpr.csv", state.observationFalsePositiveRate));
//    loggers.push_back(new ValueLogger(paramOutput / "fnr.csv", state.observationFalseNegativeRate));
//    loggers.push_back(new ValueLogger(paramOutput / "geo_gen_prob.csv", state.geometricGenerationProb));
//    loggers.push_back(new ValueLogger(paramOutput / "zt_mult_binom_prob.csv", state.ztMultiplicativeBinomialProb));
//    loggers.push_back(new ValueLogger(paramOutput / "zt_mult_binom_assoc.csv", state.ztMultiplicativeBinomialAssoc));
//    loggers.push_back(new ValueLogger(paramOutput / "geometric_coi_prob.csv", state.geometricCOIProb));
//    loggers.push_back(new ValueLogger(paramOutput / "infection_order.csv", state.infectionEventOrdering));
//    loggers.push_back(new ValueLogger(paramOutput / "likelihood.csv", model));
//    for(const auto& [locus_label, locus] : state.loci) {
//        loggers.push_back(new ValueLogger(paramOutput / (locus_label + "_frequencies.csv"), state.alleleFrequencies.alleleFrequencies(locus)));
//    }
//
//    for(const auto& logger : loggers) {
//        logger->reset();
//    }
//
//    state.observationFalsePositiveRate.saveState("state1");
//    auto currLik = model.value();
//    state.observationFalsePositiveRate.setValue(.2);
//    ASSERT_TRUE(model.isDirty());
//    auto newLik = model.value();
//    state.observationFalsePositiveRate.restoreState("state1");
//    auto oldLik = model.value();
//    std::cout << "FPR Old/New: " << oldLik << " ||| " << newLik << std::endl;
//    ASSERT_DOUBLE_EQ(currLik, oldLik);
//
//
//    state.observationFalseNegativeRate.saveState("state1");
//    currLik = model.value();
//    state.observationFalseNegativeRate.setValue(.2);
//    ASSERT_TRUE(model.isDirty());
//    newLik = model.value();
//    ASSERT_FALSE(model.isDirty());
//    state.observationFalseNegativeRate.restoreState("state1");
//    oldLik = model.value();
//    std::cout << "FNR Old/New: " << oldLik << " ||| " << newLik << std::endl;
//    ASSERT_DOUBLE_EQ(currLik, oldLik);
//
//    state.geometricGenerationProb.saveState("state1");
//    currLik = model.value();
//    state.geometricGenerationProb.setValue(.99999);
//    ASSERT_TRUE(model.isDirty());
//    newLik = model.value();
//    ASSERT_FALSE(model.isDirty());
//    state.geometricGenerationProb.restoreState("state1");
//    oldLik = model.value();
//    std::cout << "Geo Gen Prob Old/New: " << oldLik << " ||| " << newLik << std::endl;
//    ASSERT_DOUBLE_EQ(currLik, oldLik);
//
//    state.ztMultiplicativeBinomialProb.saveState("state1");
//    currLik = model.value();
//    state.ztMultiplicativeBinomialProb.setValue(.9);
//    ASSERT_TRUE(model.isDirty());
//    newLik = model.value();
//    ASSERT_FALSE(model.isDirty());
//    state.ztMultiplicativeBinomialProb.restoreState("state1");
//    oldLik = model.value();
//    std::cout << "ZTMB Prob Old/New: " << oldLik << " ||| " << newLik << std::endl;
//    ASSERT_DOUBLE_EQ(currLik, oldLik);
//
//    state.ztMultiplicativeBinomialAssoc.saveState("state1");
//    currLik = model.value();
//    state.ztMultiplicativeBinomialAssoc.setValue(.2);
//    ASSERT_TRUE(model.isDirty());
//    newLik = model.value();
//    ASSERT_FALSE(model.isDirty());
//    state.ztMultiplicativeBinomialAssoc.restoreState("state1");
//    oldLik = model.value();
//    std::cout << "ZTMB Assoc Old/New: " << oldLik << " ||| " << newLik << std::endl;
//    ASSERT_DOUBLE_EQ(currLik, oldLik);
//
//    state.geometricCOIProb.saveState("state1");
//    currLik = model.value();
//    state.geometricCOIProb.setValue(.99999);
//    ASSERT_TRUE(model.isDirty());
//    newLik = model.value();
//    ASSERT_FALSE(model.isDirty());
//    state.geometricCOIProb.restoreState("state1");
//    oldLik = model.value();
//    std::cout << "Geo COI Prob Old/New: " << oldLik << " ||| " << newLik << std::endl;
//    ASSERT_DOUBLE_EQ(currLik, oldLik);
//
//
//    boost::random::mt19937 r;
//    RandomizedScheduler scheduler(&r, 500);
//    scheduler.registerSampler(new OrderSampler(state.infectionEventOrdering, model, &r, 1));
//    scheduler.registerSampler(new OrderSampler(state.infectionEventOrdering, model, &r, 2));
//    scheduler.registerSampler(new OrderSampler(state.infectionEventOrdering, model, &r, 3));
//    scheduler.registerSampler(new OrderSampler(state.infectionEventOrdering, model, &r, 4));
//    scheduler.registerSampler(new OrderSampler(state.infectionEventOrdering, model, &r, 5));
//    scheduler.registerSampler(new OrderSampler(state.infectionEventOrdering, model, &r, 6));
//    scheduler.registerSampler(new OrderSampler(state.infectionEventOrdering, model, &r, 7));
//    scheduler.registerSampler(new OrderSampler(state.infectionEventOrdering, model, &r, 8));
//    scheduler.registerSampler(new OrderSampler(state.infectionEventOrdering, model, &r, 9));
//    scheduler.registerSampler(new ConstrainedContinuousRandomWalk(state.observationFalsePositiveRate, model, 0.0, 1.0, &r));
//    scheduler.registerSampler(new ConstrainedContinuousRandomWalk(state.observationFalseNegativeRate, model, 0.0, 1.0, &r));
//    scheduler.registerSampler(new ConstrainedContinuousRandomWalk(state.geometricGenerationProb, model, 0.0, 1.0, &r));
//    scheduler.registerSampler(new ConstrainedContinuousRandomWalk(state.ztMultiplicativeBinomialProb, model, 0.0, 1.0, &r));
//    scheduler.registerSampler(new ConstrainedContinuousRandomWalk(state.ztMultiplicativeBinomialAssoc, model, 0.0, std::numeric_limits<double>::max(), &r));
//    scheduler.registerSampler(new ConstrainedContinuousRandomWalk(state.geometricCOIProb, model, 0.0, 1.0, &r));
//
//
//    for(auto &infection : state.infections) {
//        for(const auto& [locus_label, locus] : state.loci) {
//            if (infection->latentGenotype().contains(locus)) {
//                auto &latentGenotype = infection->latentGenotype(locus);
//                scheduler.registerSampler(new GeneticsSampler(latentGenotype, model, &r));
//            }
//        }
//    }
//
//    for(const auto& [locus_label, locus] : state.loci) {
//        scheduler.registerSampler(new SALTSampler<ModelOne>(state.alleleFrequencies.alleleFrequencies(locus), model, &r));
//    }
//
//    for (int k = 0; k < 50000; ++k) {
////        scheduler.update();
//        scheduler.step();
//        if(k % 10 == 0) {
//            std::cout << "Current LLik: " << model.value() << std::endl;
//        }
//    }
//
//    for (int i = 0; i < 50000; ++i) {
//        scheduler.step();
//        if (i % 10 == 0) {
//            for (const auto& logger : loggers) {
//                logger->logValue();
//            }
//        }
//    }
//
//    std::cout << "Completed" << std::endl;
//}