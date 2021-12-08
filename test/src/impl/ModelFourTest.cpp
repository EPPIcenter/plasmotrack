////
//// Created by Maxwell Murphy on 12/10/20.
////
//
////#include <boost/filesystem.hpp>
//
//#include <filesystem>
//#include <omp.h>
//#include <memory>
//#include <random>
//
//#include "gtest/gtest.h"
//
//#include "core/io/loggers/AbstractLogger.h"
//#include "core/io/loggers/FileOutput.h"
//#include "core/io/loggers/ParentSetDistLogger.h"
//#include "core/io/loggers/ValueLogger.h"
//#include "core/io/parse_json.h"
//#include "core/io/path_parsing.h"
//#include "core/io/serialize.h"
//
//#include "core/datatypes/Simplex.h"
//
//#include "core/samplers/RandomizedScheduler.h"
//#include "core/samplers/general/ConstrainedContinuousRandomWalk.h"
//#include "core/samplers/general/SALTSampler.h"
//#include "core/samplers/genetics/RandomAllelesBitSetSampler.h"
//#include "core/samplers/genetics/ZanellaAllelesBitSetSampler.h"
//#include "core/samplers/topology/OrderSampler.h"
//#include "core/samplers/topology/ZanellaNeighborOrderSampler.h"
//#include "core/samplers/topology/ZanellaOrderSampler.h"
//
//#include "impl/model/ModelFour.h"
//
//
//namespace fs = std::filesystem;
//
//using namespace transmission_nets::impl;
//using namespace transmission_nets::core::io;
//using namespace transmission_nets::core::datatypes;
//using namespace transmission_nets::core::samplers;
//using namespace transmission_nets::core::parameters;
//
//TEST(ModelFourTest, CoreTest) {
//
//    // #pragma omp parallel
//    //     {
//    //         int ID = omp_get_thread_num();
//    //         std::cout << "Hello " << ID << std::endl;
//    //         std::cout << "World" << ID << std::endl;
//    //     }
//
//    using InfectionEvent = ModelFour::InfectionEvent;
//    using Locus = ModelFour::LocusImpl;
////    using InformedGeneticsSampler = genetics::ZanellaAllelesBitSetSampler<ModelFour::Model, boost::random::mt19937, ModelFour::GeneticsImpl>;
//    //    using GeneticsSampler = genetics::RandomAllelesBitSetSampler<ModelFour, boost::random::mt19937, ModelFour::GeneticsImpl>;
//
//    const auto testsDir = getPathFromEnvVar("TRANSMISSION_NETWORK_TESTS_DIR");
//    const fs::path nodesFile{"resources/JSON/nodes6.json"};
//    const fs::path outputDir{"outputs/ModelFourTests/CoreTest"};
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
//    ModelFour::Model model(loci, infections, disallowedParents);
//    ModelFour::StateLogger logger(model, testsDir / outputDir);
//
//    boost::random::mt19937 r;
//    RandomizedScheduler scheduler(&r, 50000);
//
//    scheduler.registerSampler({.sampler = new ConstrainedContinuousRandomWalk(model.state.geometricGenerationProb, model, 0.0, 1.0, &r, .1, .01, 2),
//                               //                               .adaptationStart = 1000,
//                               //                               .adaptationEnd = 3000,
//                               //                               .scaledAdaptation = true,
//                               .weight = 40});
//
//    scheduler.registerSampler({.sampler = new ConstrainedContinuousRandomWalk(model.state.lossProb, model, 0.0, 1.0, &r, .1, .01, 2),
//                               //                               .adaptationStart = 1000,
//                               //                               .adaptationEnd = 3000,
//                               //                               .updateStart = 100,
//                               //                               .scaledAdaptation = true,
//                               .weight = 40});
//
//        scheduler.registerSampler({.sampler = new ConstrainedContinuousRandomWalk(model.state.mutationProb, model, 0.0, .05, &r, .1, .01, 2),
////                                   .adaptationStart = 1000,
////                                   .adaptationEnd = 3000,
//                                   .weight = 40});
//
//    scheduler.registerSampler({.sampler = new ConstrainedContinuousRandomWalk(model.state.meanCOI, model, 0.0, double(ModelFour::MAX_COI), &r, .1, .01, 1),
//                               //                               .adaptationStart = 1000,
//                               //                               .adaptationEnd = 3000,
//                               //                               .scaledAdaptation = true,
//                               .weight = 40});
//
//    //    scheduler.registerSampler({.sampler = new ZanellaOrderSampler(model.state.infectionEventOrdering, model, &r, infections.size()),
//    ////                               .updateStart = 0,
//    //                               .weight = 1});
//    ////
//    //    scheduler.registerSampler({.sampler = new ZanellaNeighborOrderSampler(model.state.infectionEventOrdering, model, &r),
//    ////                               .updateStart = 0,
//    //                               .weight = 1});
//
//    scheduler.registerSampler({.sampler = new topology::OrderSampler(model.state.infectionEventOrdering, model, &r, infections.size() / 2),
//                               //                           .adaptationStart = 1000,
//                               //                           .adaptationEnd = 3000,
//                               .weight = 1000});
//
//    //    for (int l = 1; l < (int) model.state.infections.size() / 2; ++l) {
//    //        scheduler.registerSampler({.sampler = new OrderSampler(model.state.infectionEventOrdering, model, &r, 5),
//    //                                   .adaptationStart = 1000,
//    //                                   .adaptationEnd = 3000,
//    //                                   .weight = 400});
//    //    }
//
//    for (auto &infection : model.state.infections) {
//        for (const auto &[locus_label, locus] : model.state.loci) {
//            if (infection->latentGenotype().contains(locus)) {
//                auto &latentGenotype = infection->latentGenotype(locus);
//                //                scheduler.registerSampler({.sampler = new genetics::ZanellaAllelesBitSetSampler(latentGenotype, model, &r),
//                ////                                           .updateStart = 0,
//                ////                                           .updateEnd = 1000,
//                //                                                  .weight = .02});
//                scheduler.registerSampler({.sampler = new genetics::RandomAllelesBitSetSampler(latentGenotype, model, &r),
//                                           //                                           .updateStart = 0,
//                                           .weight = 2});
//            }
//        }
//    }
//
//    for (auto &infFNR : model.state.observationFalseNegativeRates) {
//        scheduler.registerSampler({.sampler = new ConstrainedContinuousRandomWalk(infFNR, model, 0.0, 0.4, &r, .1),
//                                   //                               .adaptationStart = 1000,
//                                   //                               .adaptationEnd = 3000,
//                                   //                                          .updateStart = 1000,
//                                   .weight = 6});
//    }
//
//    for (auto &infFPR : model.state.observationFalsePositiveRates) {
//        scheduler.registerSampler({.sampler = new ConstrainedContinuousRandomWalk(infFPR, model, 0.0, 0.4, &r, .1),
//                                   //                               .adaptationStart = 1000,
//                                   //                               .adaptationEnd = 3000,
//                                   //                                          .updateStart = 1000,
//                                   .weight = 6});
//    }
//
//    for (const auto &[locus_label, locus] : model.state.loci) {
//        scheduler.registerSampler({.sampler = new SALTSampler(model.state.alleleFrequencies.alleleFrequencies(locus), model, &r, .5, .01, .5),
//                                   //                                   .adaptationStart = 100,
//                                   //                                   .adaptationEnd = 3000,
////                                   .updateStart = 1000,
//                                   .weight = 1});
//    }
//
//    std::cout << "Current LLik: " << model.value() << std::endl;
//    for (int k = 0; k < 3; ++k) {
//        scheduler.step();
//        std::cout << "(k=" << k << ") Current LLik: " << model.value() << std::endl;
//        if (k % 1 == 0) {
//            std::cout << "Current LLik: " << model.value() << std::endl;
//            std::cout << serialize(model.state.infectionEventOrdering.value()) << std::endl;
//            std::cout << serialize(model.state.alleleFrequencies.alleleFrequencies(model.state.loci.begin()->second).value()) << std::endl;
//            std::cout << model.state.meanCOI.value() << " ";
//            std::cout << model.state.geometricGenerationProb.value() << " ";
//            std::cout << model.state.observationFalsePositiveRates[0].value() << " ";
//            std::cout << model.state.observationFalseNegativeRates[0].value() << " ";
//            std::cout << model.state.lossProb.value() << std::endl;
//            std::cout << model.state.mutationProb.value() << std::endl;
//            std::cout << model.transmissionProcessList[0]->child_.id() << ": ";
//            auto distr = model.transmissionProcessList[0]->calcParentSetDist();
//            for (const auto &[llik, ps] : distr.parentSetLliks) {
//                std::cout << serialize(ps) << ": " << std::exp(llik - distr.totalLlik) << std::endl;
//            }
//            std::cout << "{S}: " << std::exp(distr.sourceLlik - distr.totalLlik) << std::endl;
//        }
//    }
//
//    for (int i = 0; i < 50000; ++i) {
//        scheduler.step();
//        std::cout << "(i=" << i << ") Current LLik: " << model.value() << std::endl;
//        if (i % 1 == 0) {
//            std::cout << "Current LLik: " << model.value() << std::endl;
//            std::cout << serialize(model.state.infectionEventOrdering.value()) << std::endl;
//            std::cout << serialize(model.state.alleleFrequencies.alleleFrequencies(model.state.loci.begin()->second).value()) << std::endl;
//            std::cout << model.state.meanCOI.value() << " ";
//            std::cout << model.state.geometricGenerationProb.value() << " ";
//            std::cout << model.state.observationFalsePositiveRates[0].value() << " ";
//            std::cout << model.state.observationFalseNegativeRates[0].value() << " ";
//            std::cout << model.state.lossProb.value() << " ";
//            std::cout << model.state.mutationProb.value() << std::endl;
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
