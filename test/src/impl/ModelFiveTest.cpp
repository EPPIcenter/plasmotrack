////
//// Created by Maxwell Murphy on 4/13/21.
////
//
////#include <boost/filesystem.hpp>
//
//#include <filesystem>
//#include <random>
//
//#include "gtest/gtest.h"
//
//#include "core/io/loggers/AbstractLogger.h"
//#include "core/io/parse_json.h"
//#include "core/io/serialize.h"
//#include "core/io/utils.h"
//#include "core/samplers/meta/ReplicaExchange.h"
//#include "impl/model/ModelFive/Model.h"
//#include "impl/model/ModelFive/ModelLogger.h"
//#include "impl/model/ModelFive/SampleScheduler.h"
//#include "impl/model/ModelFive/State.h"
//#include "impl/model/ModelFive/StateLogger.h"
//
//
//namespace fs = std::filesystem;
//
//using namespace transmission_nets::impl;
//using namespace transmission_nets::core::io;
//using namespace transmission_nets::core::computation;
//using namespace transmission_nets::core::samplers;
//
//TEST(ModelFiveTest, CoreTestFull) {
//
//    const auto testsDir = getPathFromEnvVar("TRANSMISSION_NETWORK_TESTS_DIR");
////    const fs::path nodesFile{"resources/JSON/namibia/nodes.json"};
////    const fs::path outputDir{"outputs/namibia/CoreTestFull"};
//    const fs::path nodesFile{"resources/JSON/zanzibar/nodes.json"};
//    const fs::path outputDir{"outputs/zanzibar/CoreTestFull"};
//
////    const fs::path nodesFile{"resources/JSON/obs_exp/nodes_full.json"};
////    const fs::path outputDir{"outputs/obs_exp/CoreTestFull"};
//
//
//    const fs::path testFilePath = testsDir / nodesFile;
//
//    if (!fs::exists(testFilePath)) {
//        std::cerr << "Nodes test file does not exist." << std::endl;
//        exit(1);
//    }
//
//    std::ifstream testFile{testFilePath};
//
//    if (!testFile) {
//        std::cout << "Cannot open file." << std::endl;
//        exit(1);
//    }
//
//    auto j = loadJSON(testFile);
//
//    auto r = std::make_shared<boost::random::mt19937>();
//    auto repex = ReplicaExchange<ModelFive::State, ModelFive::Model, ModelFive::SampleScheduler, ModelFive::ModelLogger, ModelFive::StateLogger>(1, 1000, 1.005, r, testsDir / outputDir, false, j);
//
//    for (int k = 0; k < 10000; ++k) {
//        repex.sample();
//        repex.logModel();
//        repex.logState();
//        std::cout << "(k=" << k << ") Current LLik: " << repex.hotValue() << std::endl;
//        if (k % 1 == 0) {
////            auto idx = repex.hot_idx_;
//            auto idx = 0;
//            std::cout << serialize(repex.chains[idx].state->infectionEventOrdering->value()) << std::endl;
//            std::cout << serialize(repex.chains[idx].state->alleleFrequencies->alleleFrequencies(repex.chains[idx].state->loci.begin()->second)->value()) << std::endl;
//            std::cout
//                << repex.chains[idx].state->meanCOI->value() << " "
//                      << repex.chains[idx].state->geometricGenerationProb->value() << " "
//                      << repex.chains[idx].state->observationFalsePositiveRates[0]->value() << " "
//                      << repex.chains[idx].state->observationFalseNegativeRates[0]->value() << " "
//                      << repex.chains[idx].state->lossProb->value() << " "
//                      << repex.chains[idx].state->mutationProb->value() << " "
//                << std::endl;
//            std::cout
//                        << repex.chains[idx].state->infections[0]->id() << " "
//                      << repex.chains[idx].state->infections[0]->infectionDuration()->value() << " "
//                      << repex.chains[idx].state->infections[0]->infectionTime() << " "
//                      << repex.chains[idx].state->infections[0]->observationTime()->value() << " "
//                                  << std::endl;
//            std::cout << repex.chains[idx].model->transmissionProcessList[0]->child_->id() << ": ";
//            auto distr = repex.chains[idx].model->transmissionProcessList[0]->calcParentSetDist();
//            for (const auto &[llik, ps] : distr.parentSetLliks) {
//                std::cout << serialize(ps) << ": " << std::exp(llik - distr.totalLlik) << std::endl;
//            }
//            std::cout << "{S}: " << std::exp(distr.sourceLlik - distr.totalLlik) << std::endl;
//        }
//    }
//
//    repex.finalize();
//}
//
//
//TEST(ModelFiveTest, TestStructure) {
//    const auto testsDir = getPathFromEnvVar("TRANSMISSION_NETWORK_TESTS_DIR");
//    const fs::path nodesFile{"resources/JSON/zanzibar/nodes.json"};
//    const fs::path outputDir{"outputs/zanzibar/CoreTestFull"};
//    const fs::path testFilePath = testsDir / nodesFile;
//
//    if (!fs::exists(testFilePath)) {
//        std::cerr << "Nodes test file does not exist." << std::endl;
//        exit(1);
//    }
//
//    std::ifstream testFile{testFilePath};
//    if (!testFile) {
//        std::cout << "Cannot open file." << std::endl;
//        exit(1);
//    }
//
//    auto j = loadJSON(testFile);
//    auto r = std::make_shared<boost::random::mt19937>();
//
//    ModelFive::State state(j);
//    ModelFive::Model model(state);
//
//    auto curr_llik = model.value();
//    state.meanCOI->saveState("1");
//    state.meanCOI->setValue(1);
//    ASSERT_TRUE(model.isDirty());
//    ASSERT_TRUE(model.sourceTransmissionProcessList.back()->isDirty());
//    ASSERT_TRUE(model.likelihood.isDirty());
//    ASSERT_TRUE(model.transmissionProcessList.back()->isDirty());
//    ASSERT_FALSE(model.nodeTransmissionProcess->isDirty());
//    auto new_llik = model.value();
//    ASSERT_FALSE(model.sourceTransmissionProcessList.back()->isDirty());
//    ASSERT_FALSE(model.likelihood.isDirty());
//    ASSERT_FALSE(model.isDirty());
//    ASSERT_FALSE(model.transmissionProcessList.back()->isDirty());
//    ASSERT_NE(curr_llik, new_llik);
//    state.meanCOI->restoreState("1");
//    ASSERT_EQ(curr_llik, model.value());
//
//}

//TEST(ModelFiveTest, CoreTest90) {
//    using InfectionEvent = ModelFive::InfectionEvent;
//    using Locus = ModelFive::LocusImpl;
//
//    const auto testsDir = getPathFromEnvVar("TRANSMISSION_NETWORK_TESTS_DIR");
//    const fs::path nodesFile{"resources/JSON/locus_subset_exp/nodes_90.json"};
//    const fs::path outputDir{"outputs/locus_subset_exp/CoreTest90"};
//
//    const fs::path testFilePath = testsDir / nodesFile;
//
//    if (!fs::exists(testFilePath)) {
//        std::cerr << "Nodes test file does not exist." << std::endl;
//        exit(1);
//    }
//
//    std::ifstream testFile{testFilePath};
//
//    if (!testFile) {
//        std::cout << "Cannot open file." << std::endl;
//        exit(1);
//    }
//
//    auto j = loadJSON(testFile);
//    auto loci = parseLociFromJSON<Locus>(j);
//    auto infections = parseInfectionsFromJSON<InfectionEvent, Locus>(j, loci);
//    auto disallowedParents = parseDisallowedParentsFromJSON(j, infections);
//
//    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//    auto rng = std::default_random_engine{seed};
//    std::shuffle(infections.begin(), infections.end(), rng);
//
//    ModelFive::Model model(loci, infections, disallowedParents);
//    ModelFive::StateLogger logger(model, testsDir / outputDir);
//
//    boost::random::mt19937 r;
//    RandomizedScheduler scheduler(&r, 10000);
//
//    long double totalInfections = model.state.infections.size();
//    long double totalLoci = model.state.loci.size();
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.geometricGenerationProb, model, 0.0, 1.0, &r, .1, .01, 2),
//                                      .weight = totalInfections
//                              });
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.lossProb, model, 0.0, 1.0, &r, .1, .01, 2),
//                                      .weight = totalInfections
//                              });
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.mutationProb, model, 0.0, .05, &r, .1, .01, 2),
//                                      .weight = totalInfections
//                              });
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.meanCOI, model, 0.0, double(ModelFive::MAX_COI), &r, .1, .01, 1),
//                                      .weight = totalInfections
//                              });
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.infectionDurationShape, model, 0.0, 100.0, &r, .1),
//                                      .weight = totalInfections
//                              });
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.infectionDurationScale, model, 0.0, 10000.0, &r, .1),
//                                      .weight = totalInfections
//                              });
//
//
//    for (auto &infection : model.state.infections) {
//        scheduler.registerSampler({
//                                          .sampler = new ConstrainedContinuousRandomWalk(infection->infectionDuration(), model, 0.0, 10000.0,  &r, 1, .1, 100),
//                                          //.adaptationEnd = 2000,
//                                          .weight = totalLoci * 10
//                                  });
//        for (const auto &[locus_label, locus] : model.state.loci) {
//            if (infection->latentGenotype().contains(locus)) {
//                auto& latentGenotype = infection->latentGenotype(locus);
//                scheduler.registerSampler({
//                                                  .sampler = new genetics::RandomAllelesBitSetSampler(latentGenotype, model, &r),
//                                                  .weight = 100
//                                          });
//                scheduler.registerSampler({
//                                                  .sampler = new genetics::ZanellaAllelesBitSetSampler(latentGenotype, model, &r),
//                                                  .weight = 1
//                                          });
//            }
//        }
//    }
//
//    for (auto &infFNR : model.state.observationFalseNegativeRates) {
//        scheduler.registerSampler({
//                                          .sampler = new ConstrainedContinuousRandomWalk(infFNR, model, 0.0, 0.4, &r, .1),
//                                          .weight = totalLoci
//                                  });
//    }
//
//    for (auto &infFPR : model.state.observationFalsePositiveRates) {
//        scheduler.registerSampler({
//                                          .sampler = new ConstrainedContinuousRandomWalk(infFPR, model, 0.0, 0.4, &r, .1),
//                                          .weight = totalLoci
//                                  });
//    }
//
//    for (const auto &[locus_label, locus] : model.state.loci) {
//        scheduler.registerSampler({
//                                          .sampler = new SALTSampler(model.state.alleleFrequencies.alleleFrequencies(locus), model, &r, 1, .01, 10),
//                                          .weight = totalInfections
//                                  });
//    }
//
//    std::cout << "Current LLik: " << model.value() << std::endl;
//    for (int k = 0; k < 5000; ++k) {
//        scheduler.step();
//        std::cout << "(k=" << k << ") Current LLik: " << model.value() << std::endl;
//        if (k % 1 == 0) {
//            std::cout << "Current LLik: " << model.value() << std::endl;
//            std::cout << serialize(model.state.infectionEventOrdering.value()) << std::endl;
//            std::cout << serialize(model.state.alleleFrequencies.alleleFrequencies(model.state.loci.begin()->second).value()) << std::endl;
//            std::cout << model.state.meanCOI.value() << " "
//                      << model.state.geometricGenerationProb.value() << " "
//                      << model.state.observationFalsePositiveRates[0].value() << " "
//                      << model.state.observationFalseNegativeRates[0].value() << " "
//                      << model.state.lossProb.value() << " "
//                      << model.state.mutationProb.value() << std::endl;
//            std::cout << model.state.infections[0]->id() << " "
//                      << model.state.infections[0]->infectionDuration().value() << " "
//                      << model.state.infections[0]->infectionTime() << " "
//                      << model.state.infections[0]->observationTime().value() << std::endl;
//            std::cout << model.transmissionProcessList[0]->child_.id() << ": ";
//            auto distr = model.transmissionProcessList[0]->calcParentSetDist();
//            for (const auto &[llik, ps] : distr.parentSetLliks) {
//                std::cout << serialize(ps) << ": " << std::exp(llik - distr.totalLlik) << std::endl;
//            }
//            std::cout << "{S}: " << std::exp(distr.sourceLlik - distr.totalLlik) << std::endl;
//        }
//    }
//
//    for (int i = 0; i < 5000; ++i) {
//        scheduler.step();
//        std::cout << "(i=" << i << ") Current LLik: " << model.value() << std::endl;
//        if (i % 1 == 0) {
//            std::cout << "Current LLik: " << model.value() << std::endl;
//            std::cout << serialize(model.state.infectionEventOrdering.value()) << std::endl;
//            std::cout << serialize(model.state.alleleFrequencies.alleleFrequencies(model.state.loci.begin()->second).value()) << std::endl;
//            std::cout << model.state.meanCOI.value() << " "
//                      << model.state.geometricGenerationProb.value() << " "
//                      << model.state.observationFalsePositiveRates[0].value() << " "
//                      << model.state.observationFalseNegativeRates[0].value() << " "
//                      << model.state.lossProb.value() << " "
//                      << model.state.mutationProb.value() << std::endl;
//            std::cout << model.state.infections[0]->id() << " "
//                      << model.state.infections[0]->infectionDuration().value() << " "
//                      << model.state.infections[0]->infectionTime() << " "
//                      << model.state.infections[0]->observationTime().value() << std::endl;
//            std::cout << model.transmissionProcessList[0]->child_.id() << ": ";
//            auto distr = model.transmissionProcessList[0]->calcParentSetDist();
//            for (const auto &[llik, ps] : distr.parentSetLliks) {
//                std::cout << serialize(ps) << ": " << std::exp(llik - distr.totalLlik) << std::endl;
//            }
//            std::cout << "{S}: " << std::exp(distr.sourceLlik - distr.totalLlik) << std::endl;
//            logger.logParameters();
//        }
//    }
//
//    std::cout << "Completed" << std::endl;
//}
//
//TEST(ModelFiveTest, CoreTest75) {
//    using InfectionEvent = ModelFive::InfectionEvent;
//    using Locus = ModelFive::LocusImpl;
//
//    const auto testsDir = getPathFromEnvVar("TRANSMISSION_NETWORK_TESTS_DIR");
//    const fs::path nodesFile{"resources/JSON/locus_subset_exp/nodes_75.json"};
//    const fs::path outputDir{"outputs/locus_subset_exp/CoreTest75"};
//
//    const fs::path testFilePath = testsDir / nodesFile;
//
//    if (!fs::exists(testFilePath)) {
//        std::cerr << "Nodes test file does not exist." << std::endl;
//        exit(1);
//    }
//
//    std::ifstream testFile{testFilePath};
//
//    if (!testFile) {
//        std::cout << "Cannot open file." << std::endl;
//        exit(1);
//    }
//
//    auto j = loadJSON(testFile);
//    auto loci = parseLociFromJSON<Locus>(j);
//    auto infections = parseInfectionsFromJSON<InfectionEvent, Locus>(j, loci);
//    auto disallowedParents = parseDisallowedParentsFromJSON(j, infections);
//
//    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//    auto rng = std::default_random_engine{seed};
//    std::shuffle(infections.begin(), infections.end(), rng);
//
//    ModelFive::Model model(loci, infections, disallowedParents);
//    ModelFive::StateLogger logger(model, testsDir / outputDir);
//
//    boost::random::mt19937 r;
//    RandomizedScheduler scheduler(&r, 10000);
//
//    long double totalInfections = model.state.infections.size();
//    long double totalLoci = model.state.loci.size();
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.geometricGenerationProb, model, 0.0, 1.0, &r, .1, .01, 2),
//                                      .weight = totalInfections
//                              });
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.lossProb, model, 0.0, 1.0, &r, .1, .01, 2),
//                                      .weight = totalInfections
//                              });
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.mutationProb, model, 0.0, .05, &r, .1, .01, 2),
//                                      .weight = totalInfections
//                              });
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.meanCOI, model, 0.0, double(ModelFive::MAX_COI), &r, .1, .01, 1),
//                                      .weight = totalInfections
//                              });
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.infectionDurationShape, model, 0.0, 100.0, &r, .1),
//                                      .weight = totalInfections
//                              });
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.infectionDurationScale, model, 0.0, 10000.0, &r, .1),
//                                      .weight = totalInfections
//                              });
//
//
//    for (auto &infection : model.state.infections) {
//        scheduler.registerSampler({
//                                          .sampler = new ConstrainedContinuousRandomWalk(infection->infectionDuration(), model, 0.0, 10000.0,  &r, 1, .1, 100),
//                                          //.adaptationEnd = 2000,
//                                          .weight = totalLoci * 10
//                                  });
//        for (const auto &[locus_label, locus] : model.state.loci) {
//            if (infection->latentGenotype().contains(locus)) {
//                auto& latentGenotype = infection->latentGenotype(locus);
//                scheduler.registerSampler({
//                                                  .sampler = new genetics::RandomAllelesBitSetSampler(latentGenotype, model, &r),
//                                                  .weight = 100
//                                          });
//                scheduler.registerSampler({
//                                                  .sampler = new genetics::ZanellaAllelesBitSetSampler(latentGenotype, model, &r),
//                                                  .weight = 1
//                                          });
//            }
//        }
//    }
//
//    for (auto &infFNR : model.state.observationFalseNegativeRates) {
//        scheduler.registerSampler({
//                                          .sampler = new ConstrainedContinuousRandomWalk(infFNR, model, 0.0, 0.4, &r, .1),
//                                          .weight = totalLoci
//                                  });
//    }
//
//    for (auto &infFPR : model.state.observationFalsePositiveRates) {
//        scheduler.registerSampler({
//                                          .sampler = new ConstrainedContinuousRandomWalk(infFPR, model, 0.0, 0.4, &r, .1),
//                                          .weight = totalLoci
//                                  });
//    }
//
//    for (const auto &[locus_label, locus] : model.state.loci) {
//        scheduler.registerSampler({
//                                          .sampler = new SALTSampler(model.state.alleleFrequencies.alleleFrequencies(locus), model, &r, 1, .01, 10),
//                                          .weight = totalInfections
//                                  });
//    }
//
//    std::cout << "Current LLik: " << model.value() << std::endl;
//    for (int k = 0; k < 5000; ++k) {
//        scheduler.step();
//        std::cout << "(k=" << k << ") Current LLik: " << model.value() << std::endl;
//        if (k % 1 == 0) {
//            std::cout << "Current LLik: " << model.value() << std::endl;
//            std::cout << serialize(model.state.infectionEventOrdering.value()) << std::endl;
//            std::cout << serialize(model.state.alleleFrequencies.alleleFrequencies(model.state.loci.begin()->second).value()) << std::endl;
//            std::cout << model.state.meanCOI.value() << " "
//                      << model.state.geometricGenerationProb.value() << " "
//                      << model.state.observationFalsePositiveRates[0].value() << " "
//                      << model.state.observationFalseNegativeRates[0].value() << " "
//                      << model.state.lossProb.value() << " "
//                      << model.state.mutationProb.value() << std::endl;
//            std::cout << model.state.infections[0]->id() << " "
//                      << model.state.infections[0]->infectionDuration().value() << " "
//                      << model.state.infections[0]->infectionTime() << " "
//                      << model.state.infections[0]->observationTime().value() << std::endl;
//            std::cout << model.transmissionProcessList[0]->child_.id() << ": ";
//            auto distr = model.transmissionProcessList[0]->calcParentSetDist();
//            for (const auto &[llik, ps] : distr.parentSetLliks) {
//                std::cout << serialize(ps) << ": " << std::exp(llik - distr.totalLlik) << std::endl;
//            }
//            std::cout << "{S}: " << std::exp(distr.sourceLlik - distr.totalLlik) << std::endl;
//        }
//    }
//
//    for (int i = 0; i < 5000; ++i) {
//        scheduler.step();
//        std::cout << "(i=" << i << ") Current LLik: " << model.value() << std::endl;
//        if (i % 1 == 0) {
//            std::cout << "Current LLik: " << model.value() << std::endl;
//            std::cout << serialize(model.state.infectionEventOrdering.value()) << std::endl;
//            std::cout << serialize(model.state.alleleFrequencies.alleleFrequencies(model.state.loci.begin()->second).value()) << std::endl;
//            std::cout << model.state.meanCOI.value() << " "
//                      << model.state.geometricGenerationProb.value() << " "
//                      << model.state.observationFalsePositiveRates[0].value() << " "
//                      << model.state.observationFalseNegativeRates[0].value() << " "
//                      << model.state.lossProb.value() << " "
//                      << model.state.mutationProb.value() << std::endl;
//            std::cout << model.state.infections[0]->id() << " "
//                      << model.state.infections[0]->infectionDuration().value() << " "
//                      << model.state.infections[0]->infectionTime() << " "
//                      << model.state.infections[0]->observationTime().value() << std::endl;
//            std::cout << model.transmissionProcessList[0]->child_.id() << ": ";
//            auto distr = model.transmissionProcessList[0]->calcParentSetDist();
//            for (const auto &[llik, ps] : distr.parentSetLliks) {
//                std::cout << serialize(ps) << ": " << std::exp(llik - distr.totalLlik) << std::endl;
//            }
//            std::cout << "{S}: " << std::exp(distr.sourceLlik - distr.totalLlik) << std::endl;
//            logger.logParameters();
//        }
//    }
//
//    std::cout << "Completed" << std::endl;
//}

//TEST(ModelFiveTest, CoreTest50) {
//    using InfectionEvent = ModelFive::InfectionEvent;
//    using Locus = ModelFive::LocusImpl;
//
//    const auto testsDir = getPathFromEnvVar("TRANSMISSION_NETWORK_TESTS_DIR");
//    const fs::path nodesFile{"resources/JSON/locus_subset_exp/nodes_50.json"};
//    const fs::path outputDir{"outputs/locus_subset_exp/CoreTest50"};
//
//    const fs::path testFilePath = testsDir / nodesFile;
//
//    if (!fs::exists(testFilePath)) {
//        std::cerr << "Nodes test file does not exist." << std::endl;
//        exit(1);
//    }
//
//    std::ifstream testFile{testFilePath};
//
//    if (!testFile) {
//        std::cout << "Cannot open file." << std::endl;
//        exit(1);
//    }
//
//    auto j = loadJSON(testFile);
//    auto loci = parseLociFromJSON<Locus>(j);
//    auto infections = parseInfectionsFromJSON<InfectionEvent, Locus>(j, loci);
//    auto disallowedParents = parseDisallowedParentsFromJSON(j, infections);
//
//    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//    auto rng = std::default_random_engine{seed};
//    std::shuffle(infections.begin(), infections.end(), rng);
//
//    ModelFive::Model model(loci, infections, disallowedParents);
//    ModelFive::StateLogger logger(model, testsDir / outputDir);
//
//    boost::random::mt19937 r;
//    RandomizedScheduler scheduler(&r, 10000);
//
//    long double totalInfections = model.state.infections.size();
//    long double totalLoci = model.state.loci.size();
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.geometricGenerationProb, model, 0.0, 1.0, &r, .1, .01, 2),
//                                      .weight = totalInfections
//                              });
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.lossProb, model, 0.0, 1.0, &r, .1, .01, 2),
//                                      .weight = totalInfections
//                              });
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.mutationProb, model, 0.0, .05, &r, .1, .01, 2),
//                                      .weight = totalInfections
//                              });
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.meanCOI, model, 0.0, double(ModelFive::MAX_COI), &r, .1, .01, 1),
//                                      .weight = totalInfections
//                              });
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.infectionDurationShape, model, 0.0, 100.0, &r, .1),
//                                      .weight = totalInfections
//                              });
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.infectionDurationScale, model, 0.0, 10000.0, &r, .1),
//                                      .weight = totalInfections
//                              });
//
//
//    for (auto &infection : model.state.infections) {
//        scheduler.registerSampler({
//                                          .sampler = new ConstrainedContinuousRandomWalk(infection->infectionDuration(), model, 0.0, 10000.0,  &r, 1, .1, 100),
//                                          //.adaptationEnd = 2000,
//                                          .weight = totalLoci * 10
//                                  });
//        for (const auto &[locus_label, locus] : model.state.loci) {
//            if (infection->latentGenotype().contains(locus)) {
//                auto& latentGenotype = infection->latentGenotype(locus);
//                scheduler.registerSampler({
//                                                  .sampler = new genetics::RandomAllelesBitSetSampler(latentGenotype, model, &r),
//                                                  .weight = 100
//                                          });
//                scheduler.registerSampler({
//                                                  .sampler = new genetics::ZanellaAllelesBitSetSampler(latentGenotype, model, &r),
//                                                  .weight = 1
//                                          });
//            }
//        }
//    }
//
//    for (auto &infFNR : model.state.observationFalseNegativeRates) {
//        scheduler.registerSampler({
//                                          .sampler = new ConstrainedContinuousRandomWalk(infFNR, model, 0.0, 0.4, &r, .1),
//                                          .weight = totalLoci
//                                  });
//    }
//
//    for (auto &infFPR : model.state.observationFalsePositiveRates) {
//        scheduler.registerSampler({
//                                          .sampler = new ConstrainedContinuousRandomWalk(infFPR, model, 0.0, 0.4, &r, .1),
//                                          .weight = totalLoci
//                                  });
//    }
//
//    for (const auto &[locus_label, locus] : model.state.loci) {
//        scheduler.registerSampler({
//                                          .sampler = new SALTSampler(model.state.alleleFrequencies.alleleFrequencies(locus), model, &r, 1, .01, 10),
//                                          .weight = totalInfections
//                                  });
//    }
//
//    std::cout << "Current LLik: " << model.value() << std::endl;
//    for (int k = 0; k < 5000; ++k) {
//        scheduler.step();
//        std::cout << "(k=" << k << ") Current LLik: " << model.value() << std::endl;
//        if (k % 1 == 0) {
//            std::cout << "Current LLik: " << model.value() << std::endl;
//            std::cout << serialize(model.state.infectionEventOrdering.value()) << std::endl;
//            std::cout << serialize(model.state.alleleFrequencies.alleleFrequencies(model.state.loci.begin()->second).value()) << std::endl;
//            std::cout << model.state.meanCOI.value() << " "
//                      << model.state.geometricGenerationProb.value() << " "
//                      << model.state.observationFalsePositiveRates[0].value() << " "
//                      << model.state.observationFalseNegativeRates[0].value() << " "
//                      << model.state.lossProb.value() << " "
//                      << model.state.mutationProb.value() << std::endl;
//            std::cout << model.state.infections[0]->id() << " "
//                      << model.state.infections[0]->infectionDuration().value() << " "
//                      << model.state.infections[0]->infectionTime() << " "
//                      << model.state.infections[0]->observationTime().value() << std::endl;
//            std::cout << model.transmissionProcessList[0]->child_.id() << ": ";
//            auto distr = model.transmissionProcessList[0]->calcParentSetDist();
//            for (const auto &[llik, ps] : distr.parentSetLliks) {
//                std::cout << serialize(ps) << ": " << std::exp(llik - distr.totalLlik) << std::endl;
//            }
//            std::cout << "{S}: " << std::exp(distr.sourceLlik - distr.totalLlik) << std::endl;
//        }
//    }
//
//    for (int i = 0; i < 5000; ++i) {
//        scheduler.step();
//        std::cout << "(i=" << i << ") Current LLik: " << model.value() << std::endl;
//        if (i % 1 == 0) {
//            std::cout << "Current LLik: " << model.value() << std::endl;
//            std::cout << serialize(model.state.infectionEventOrdering.value()) << std::endl;
//            std::cout << serialize(model.state.alleleFrequencies.alleleFrequencies(model.state.loci.begin()->second).value()) << std::endl;
//            std::cout << model.state.meanCOI.value() << " "
//                      << model.state.geometricGenerationProb.value() << " "
//                      << model.state.observationFalsePositiveRates[0].value() << " "
//                      << model.state.observationFalseNegativeRates[0].value() << " "
//                      << model.state.lossProb.value() << " "
//                      << model.state.mutationProb.value() << std::endl;
//            std::cout << model.state.infections[0]->id() << " "
//                      << model.state.infections[0]->infectionDuration().value() << " "
//                      << model.state.infections[0]->infectionTime() << " "
//                      << model.state.infections[0]->observationTime().value() << std::endl;
//            std::cout << model.transmissionProcessList[0]->child_.id() << ": ";
//            auto distr = model.transmissionProcessList[0]->calcParentSetDist();
//            for (const auto &[llik, ps] : distr.parentSetLliks) {
//                std::cout << serialize(ps) << ": " << std::exp(llik - distr.totalLlik) << std::endl;
//            }
//            std::cout << "{S}: " << std::exp(distr.sourceLlik - distr.totalLlik) << std::endl;
//            logger.logParameters();
//        }
//    }
//
//    std::cout << "Completed" << std::endl;
//}
//
//TEST(ModelFiveTest, CoreTest25) {
//    using InfectionEvent = ModelFive::InfectionEvent;
//    using Locus = ModelFive::LocusImpl;
//
//    const auto testsDir = getPathFromEnvVar("TRANSMISSION_NETWORK_TESTS_DIR");
//    const fs::path nodesFile{"resources/JSON/locus_subset_exp/nodes_25.json"};
//    const fs::path outputDir{"outputs/locus_subset_exp/CoreTest25"};
//
//    const fs::path testFilePath = testsDir / nodesFile;
//
//    if (!fs::exists(testFilePath)) {
//        std::cerr << "Nodes test file does not exist." << std::endl;
//        exit(1);
//    }
//
//    std::ifstream testFile{testFilePath};
//
//    if (!testFile) {
//        std::cout << "Cannot open file." << std::endl;
//        exit(1);
//    }
//
//    auto j = loadJSON(testFile);
//    auto loci = parseLociFromJSON<Locus>(j);
//    auto infections = parseInfectionsFromJSON<InfectionEvent, Locus>(j, loci);
//    auto disallowedParents = parseDisallowedParentsFromJSON(j, infections);
//
//    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//    auto rng = std::default_random_engine{seed};
//    std::shuffle(infections.begin(), infections.end(), rng);
//
//    ModelFive::Model model(loci, infections, disallowedParents);
//    ModelFive::StateLogger logger(model, testsDir / outputDir);
//
//    boost::random::mt19937 r;
//    RandomizedScheduler scheduler(&r, 10000);
//
//    long double totalInfections = model.state.infections.size();
//    long double totalLoci = model.state.loci.size();
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.geometricGenerationProb, model, 0.0, 1.0, &r, .1, .01, 2),
//                                      .weight = totalInfections
//                              });
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.lossProb, model, 0.0, 1.0, &r, .1, .01, 2),
//                                      .weight = totalInfections
//                              });
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.mutationProb, model, 0.0, .05, &r, .1, .01, 2),
//                                      .weight = totalInfections
//                              });
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.meanCOI, model, 0.0, double(ModelFive::MAX_COI), &r, .1, .01, 1),
//                                      .weight = totalInfections
//                              });
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.infectionDurationShape, model, 0.0, 100.0, &r, .1),
//                                      .weight = totalInfections
//                              });
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.infectionDurationScale, model, 0.0, 10000.0, &r, .1),
//                                      .weight = totalInfections
//                              });
//
//
//    for (auto &infection : model.state.infections) {
//        scheduler.registerSampler({
//                                          .sampler = new ConstrainedContinuousRandomWalk(infection->infectionDuration(), model, 0.0, 10000.0,  &r, 1, .1, 100),
//                                          //.adaptationEnd = 2000,
//                                          .weight = totalLoci * 10
//                                  });
//        for (const auto &[locus_label, locus] : model.state.loci) {
//            if (infection->latentGenotype().contains(locus)) {
//                auto& latentGenotype = infection->latentGenotype(locus);
//                scheduler.registerSampler({
//                                                  .sampler = new genetics::RandomAllelesBitSetSampler(latentGenotype, model, &r),
//                                                  .weight = 100
//                                          });
//                scheduler.registerSampler({
//                                                  .sampler = new genetics::ZanellaAllelesBitSetSampler(latentGenotype, model, &r),
//                                                  .weight = 1
//                                          });
//            }
//        }
//    }
//
//    for (auto &infFNR : model.state.observationFalseNegativeRates) {
//        scheduler.registerSampler({
//                                          .sampler = new ConstrainedContinuousRandomWalk(infFNR, model, 0.0, 0.4, &r, .1),
//                                          .weight = totalLoci
//                                  });
//    }
//
//    for (auto &infFPR : model.state.observationFalsePositiveRates) {
//        scheduler.registerSampler({
//                                          .sampler = new ConstrainedContinuousRandomWalk(infFPR, model, 0.0, 0.4, &r, .1),
//                                          .weight = totalLoci
//                                  });
//    }
//
//    for (const auto &[locus_label, locus] : model.state.loci) {
//        scheduler.registerSampler({
//                                          .sampler = new SALTSampler(model.state.alleleFrequencies.alleleFrequencies(locus), model, &r, 1, .01, 10),
//                                          .weight = totalInfections
//                                  });
//    }
//
//    std::cout << "Current LLik: " << model.value() << std::endl;
//    for (int k = 0; k < 5000; ++k) {
//        scheduler.step();
//        std::cout << "(k=" << k << ") Current LLik: " << model.value() << std::endl;
//        if (k % 1 == 0) {
//            std::cout << "Current LLik: " << model.value() << std::endl;
//            std::cout << serialize(model.state.infectionEventOrdering.value()) << std::endl;
//            std::cout << serialize(model.state.alleleFrequencies.alleleFrequencies(model.state.loci.begin()->second).value()) << std::endl;
//            std::cout << model.state.meanCOI.value() << " "
//                      << model.state.geometricGenerationProb.value() << " "
//                      << model.state.observationFalsePositiveRates[0].value() << " "
//                      << model.state.observationFalseNegativeRates[0].value() << " "
//                      << model.state.lossProb.value() << " "
//                      << model.state.mutationProb.value() << std::endl;
//            std::cout << model.state.infections[0]->id() << " "
//                      << model.state.infections[0]->infectionDuration().value() << " "
//                      << model.state.infections[0]->infectionTime() << " "
//                      << model.state.infections[0]->observationTime().value() << std::endl;
//            std::cout << model.transmissionProcessList[0]->child_.id() << ": ";
//            auto distr = model.transmissionProcessList[0]->calcParentSetDist();
//            for (const auto &[llik, ps] : distr.parentSetLliks) {
//                std::cout << serialize(ps) << ": " << std::exp(llik - distr.totalLlik) << std::endl;
//            }
//            std::cout << "{S}: " << std::exp(distr.sourceLlik - distr.totalLlik) << std::endl;
//        }
//    }
//
//    for (int i = 0; i < 5000; ++i) {
//        scheduler.step();
//        std::cout << "(i=" << i << ") Current LLik: " << model.value() << std::endl;
//        if (i % 1 == 0) {
//            std::cout << "Current LLik: " << model.value() << std::endl;
//            std::cout << serialize(model.state.infectionEventOrdering.value()) << std::endl;
//            std::cout << serialize(model.state.alleleFrequencies.alleleFrequencies(model.state.loci.begin()->second).value()) << std::endl;
//            std::cout << model.state.meanCOI.value() << " "
//                      << model.state.geometricGenerationProb.value() << " "
//                      << model.state.observationFalsePositiveRates[0].value() << " "
//                      << model.state.observationFalseNegativeRates[0].value() << " "
//                      << model.state.lossProb.value() << " "
//                      << model.state.mutationProb.value() << std::endl;
//            std::cout << model.state.infections[0]->id() << " "
//                      << model.state.infections[0]->infectionDuration().value() << " "
//                      << model.state.infections[0]->infectionTime() << " "
//                      << model.state.infections[0]->observationTime().value() << std::endl;
//            std::cout << model.transmissionProcessList[0]->child_.id() << ": ";
//            auto distr = model.transmissionProcessList[0]->calcParentSetDist();
//            for (const auto &[llik, ps] : distr.parentSetLliks) {
//                std::cout << serialize(ps) << ": " << std::exp(llik - distr.totalLlik) << std::endl;
//            }
//            std::cout << "{S}: " << std::exp(distr.sourceLlik - distr.totalLlik) << std::endl;
//            logger.logParameters();
//        }
//    }
//
//    std::cout << "Completed" << std::endl;
//}
//
//TEST(ModelFiveTest, CoreTest125) {
//    using InfectionEvent = ModelFive::InfectionEvent;
//    using Locus = ModelFive::LocusImpl;
//
//    const auto testsDir = getPathFromEnvVar("TRANSMISSION_NETWORK_TESTS_DIR");
//    const fs::path nodesFile{"resources/JSON/locus_subset_exp/nodes_125.json"};
//    const fs::path outputDir{"outputs/locus_subset_exp/CoreTest125"};
//
//    const fs::path testFilePath = testsDir / nodesFile;
//
//    if (!fs::exists(testFilePath)) {
//        std::cerr << "Nodes test file does not exist." << std::endl;
//        exit(1);
//    }
//
//    std::ifstream testFile{testFilePath};
//
//    if (!testFile) {
//        std::cout << "Cannot open file." << std::endl;
//        exit(1);
//    }
//
//    auto j = loadJSON(testFile);
//    auto loci = parseLociFromJSON<Locus>(j);
//    auto infections = parseInfectionsFromJSON<InfectionEvent, Locus>(j, loci);
//    auto disallowedParents = parseDisallowedParentsFromJSON(j, infections);
//
//    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//    auto rng = std::default_random_engine{seed};
//    std::shuffle(infections.begin(), infections.end(), rng);
//
//    ModelFive::Model model(loci, infections, disallowedParents);
//    ModelFive::StateLogger logger(model, testsDir / outputDir);
//
//    boost::random::mt19937 r;
//    RandomizedScheduler scheduler(&r, 10000);
//
//    long double totalInfections = model.state.infections.size();
//    long double totalLoci = model.state.loci.size();
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.geometricGenerationProb, model, 0.0, 1.0, &r, .1, .01, 2),
//                                      .weight = totalInfections
//                              });
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.lossProb, model, 0.0, 1.0, &r, .1, .01, 2),
//                                      .weight = totalInfections
//                              });
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.mutationProb, model, 0.0, .05, &r, .1, .01, 2),
//                                      .weight = totalInfections
//                              });
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.meanCOI, model, 0.0, double(ModelFive::MAX_COI), &r, .1, .01, 1),
//                                      .weight = totalInfections
//                              });
//
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.infectionDurationShape, model, 0.0, 100.0, &r, .1),
//                                      .weight = totalInfections
//                              });
//    scheduler.registerSampler({
//                                      .sampler = new ConstrainedContinuousRandomWalk(model.state.infectionDurationScale, model, 0.0, 10000.0, &r, .1),
//                                      .weight = totalInfections
//                              });
//
//
//    for (auto &infection : model.state.infections) {
//        scheduler.registerSampler({
//                                          .sampler = new ConstrainedContinuousRandomWalk(infection->infectionDuration(), model, 0.0, 10000.0,  &r, 1, .1, 100),
//                                          //.adaptationEnd = 2000,
//                                          .weight = totalLoci * 10
//                                  });
//        for (const auto &[locus_label, locus] : model.state.loci) {
//            if (infection->latentGenotype().contains(locus)) {
//                auto& latentGenotype = infection->latentGenotype(locus);
//                scheduler.registerSampler({
//                                                  .sampler = new genetics::RandomAllelesBitSetSampler(latentGenotype, model, &r),
//                                                  .weight = 100
//                                          });
//                scheduler.registerSampler({
//                                                  .sampler = new genetics::ZanellaAllelesBitSetSampler(latentGenotype, model, &r),
//                                                  .weight = 1
//                                          });
//            }
//        }
//    }
//
//    for (auto &infFNR : model.state.observationFalseNegativeRates) {
//        scheduler.registerSampler({
//                                          .sampler = new ConstrainedContinuousRandomWalk(infFNR, model, 0.0, 0.4, &r, .1),
//                                          .weight = totalLoci
//                                  });
//    }
//
//    for (auto &infFPR : model.state.observationFalsePositiveRates) {
//        scheduler.registerSampler({
//                                          .sampler = new ConstrainedContinuousRandomWalk(infFPR, model, 0.0, 0.4, &r, .1),
//                                          .weight = totalLoci
//                                  });
//    }
//
//    for (const auto &[locus_label, locus] : model.state.loci) {
//        scheduler.registerSampler({
//                                          .sampler = new SALTSampler(model.state.alleleFrequencies.alleleFrequencies(locus), model, &r, 1, .01, 10),
//                                          .weight = totalInfections
//                                  });
//    }
//
//    std::cout << "Current LLik: " << model.value() << std::endl;
//    for (int k = 0; k < 5000; ++k) {
//        scheduler.step();
//        std::cout << "(k=" << k << ") Current LLik: " << model.value() << std::endl;
//        if (k % 1 == 0) {
//            std::cout << "Current LLik: " << model.value() << std::endl;
//            std::cout << serialize(model.state.infectionEventOrdering.value()) << std::endl;
//            std::cout << serialize(model.state.alleleFrequencies.alleleFrequencies(model.state.loci.begin()->second).value()) << std::endl;
//            std::cout << model.state.meanCOI.value() << " "
//                      << model.state.geometricGenerationProb.value() << " "
//                      << model.state.observationFalsePositiveRates[0].value() << " "
//                      << model.state.observationFalseNegativeRates[0].value() << " "
//                      << model.state.lossProb.value() << " "
//                      << model.state.mutationProb.value() << std::endl;
//            std::cout << model.state.infections[0]->id() << " "
//                      << model.state.infections[0]->infectionDuration().value() << " "
//                      << model.state.infections[0]->infectionTime() << " "
//                      << model.state.infections[0]->observationTime().value() << std::endl;
//            std::cout << model.transmissionProcessList[0]->child_.id() << ": ";
//            auto distr = model.transmissionProcessList[0]->calcParentSetDist();
//            for (const auto &[llik, ps] : distr.parentSetLliks) {
//                std::cout << serialize(ps) << ": " << std::exp(llik - distr.totalLlik) << std::endl;
//            }
//            std::cout << "{S}: " << std::exp(distr.sourceLlik - distr.totalLlik) << std::endl;
//        }
//    }
//
//    for (int i = 0; i < 5000; ++i) {
//        scheduler.step();
//        std::cout << "(i=" << i << ") Current LLik: " << model.value() << std::endl;
//        if (i % 1 == 0) {
//            std::cout << "Current LLik: " << model.value() << std::endl;
//            std::cout << serialize(model.state.infectionEventOrdering.value()) << std::endl;
//            std::cout << serialize(model.state.alleleFrequencies.alleleFrequencies(model.state.loci.begin()->second).value()) << std::endl;
//            std::cout << model.state.meanCOI.value() << " "
//                      << model.state.geometricGenerationProb.value() << " "
//                      << model.state.observationFalsePositiveRates[0].value() << " "
//                      << model.state.observationFalseNegativeRates[0].value() << " "
//                      << model.state.lossProb.value() << " "
//                      << model.state.mutationProb.value() << std::endl;
//            std::cout << model.state.infections[0]->id() << " "
//                      << model.state.infections[0]->infectionDuration().value() << " "
//                      << model.state.infections[0]->infectionTime() << " "
//                      << model.state.infections[0]->observationTime().value() << std::endl;
//            std::cout << model.transmissionProcessList[0]->child_.id() << ": ";
//            auto distr = model.transmissionProcessList[0]->calcParentSetDist();
//            for (const auto &[llik, ps] : distr.parentSetLliks) {
//                std::cout << serialize(ps) << ": " << std::exp(llik - distr.totalLlik) << std::endl;
//            }
//            std::cout << "{S}: " << std::exp(distr.sourceLlik - distr.totalLlik) << std::endl;
//            logger.logParameters();
//        }
//    }
//
//    std::cout << "Completed" << std::endl;
//}
