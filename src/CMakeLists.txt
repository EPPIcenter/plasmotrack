cmake_minimum_required(VERSION 3.17)
project(transmission_networks)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")

find_package(Boost 1.74 REQUIRED container math_c99l random program_options)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(nlohmann_json 3.2.0 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(fmt REQUIRED)
find_package(MPFR REQUIRED)
find_package(gmp REQUIRED)

if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS " ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} ")
endif ()


include_directories(SYSTEM ${Boost_INCLUDE_DIR})

set(SOURCE_FILES
        core/computation/PartialLikelihood.cpp
        core/computation/LogLikelihood.cpp
        core/containers/Locus.cpp
        core/utils/generators/CombinationIndicesGenerator.cpp
        model/observation_process/AlleleCounts.cpp
        core/samplers/general/ConstrainedDiscreteRandomWalk.h
        core/samplers/genetics/RandomAllelesBitSetSampler.h
        core/utils/generators/RandomSequence.h
        core/datatypes/Simplex.cpp
        core/samplers/scheduler/Scheduler.cpp
        core/samplers/scheduler/Scheduler.h
        impl/state/ModelOneState.cpp
        impl/state/ModelOneState.h
        impl/model/ModelOne.cpp
        impl/model/ModelOne.h
        core/computation/transformers/Transformer.h
        core/computation/transformers/LogTransformer.h
        core/utils/generators/CombinationsWithRepetitionsGenerator.cpp
        core/utils/generators/CombinationsWithRepetitionsGenerator.h
        core/priors/Prior.h
        core/io/parse_json.h
        core/io/loggers/FileOutput.h
        core/io/loggers/ValueLogger.h
        core/io/serialize.h
        core/io/loggers/FileOutput.cpp
        core/samplers/scheduler/RandomizedScheduler.h
        core/io/utils.cpp
        core/io/utils.h
        model/transmission_process/node_transmission_process/NoSuperInfectionMutation.h
        model/transmission_process/NetworkBasedTransmissionProcess.h
        impl/model/ModelTwo.cpp
        impl/model/ModelTwo.h
        core/distributions/ZTPoisson.h
        core/samplers/topology/RandomAddEdgeSampler.h
        core/samplers/topology/RandomRemoveEdgeSampler.h
        core/samplers/topology/RandomSwapEdgeSampler.h
        impl/model/ModelThree.cpp
        impl/model/ModelThree.h
        impl/state/ModelThreeState.cpp
        impl/state/ModelThreeState.h
        core/io/loggers/LambdaLogger.h
        core/samplers/topology/RandomReverseEdgeSampler.h
        core/utils/ProbAnyMissing.h
        model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h
        core/distributions/pdfs/BetaLogPDF.h
        core/distributions/pdfs/BetaLogPDF.cpp
        core/distributions/pdfs/GammaLogPDF.cpp
        core/distributions/pdfs/GammaLogPDF.h
        core/mcmc/Chain.h
        core/utils/ProbAnyMissing.cpp
        core/samplers/genetics/ZanellaAllelesBitSetSampler.h
        core/samplers/topology/ZanellaNeighborOrderSampler.h
        core/samplers/specialized/ZanellaJointGeneticsGraphSampler.h
        core/samplers/topology/ZanellaOrderSampler.h
        core/samplers/general/SimplexSampler.h
        core/computation/ConstrainedOrderDerivedParentSet.h
        core/io/loggers/ParentSetDistLogger.h
        impl/model/ModelFour.cpp
        impl/model/ModelFour.h
        core/io/loggers/AbstractLogger.h
        core/io/loggers/AbstractOutput.h
        model/OrderingLikelihood.h
        core/computation/ObservationTimeDerivedOrdering.h
        impl/model/ModelFive.cpp
        impl/model/ModelFive.h
        core/utils/LogPQ.cpp
        core/utils/LogPQ.h
        core/io/loggers/MultiValueLogger.h
        core/computation/transformers/Tempered.h
        core/samplers/meta/ReplicaExchange.h

        impl/model/ModelFive/State.cpp
        impl/model/ModelFive/State.h
        impl/model/ModelFive/Model.cpp
        impl/model/ModelFive/Model.h
        impl/model/ModelFive/config.h
        impl/model/ModelFive/SampleScheduler.h
        impl/model/ModelFive/StateLogger.cpp
        impl/model/ModelFive/StateLogger.h
        impl/model/ModelFive/ModelLogger.cpp
        impl/model/ModelFive/ModelLogger.h
        impl/model/ModelFive/SamplerLogger.cpp
        impl/model/ModelFive/SamplerLogger.h

        impl/model/ModelSix/State.cpp
        impl/model/ModelSix/State.h
        impl/model/ModelSix/Model.cpp
        impl/model/ModelSix/Model.h
        impl/model/ModelSix/config.h
        impl/model/ModelSix/SampleScheduler.h
        impl/model/ModelSix/StateLogger.cpp
        impl/model/ModelSix/StateLogger.h
        impl/model/ModelSix/ModelLogger.cpp
        impl/model/ModelSix/ModelLogger.h
        impl/model/ModelSix/SamplerLogger.cpp
        impl/model/ModelSix/SamplerLogger.h

        impl/model/ModelSeven/State.cpp
        impl/model/ModelSeven/State.h
        impl/model/ModelSeven/Model.cpp
        impl/model/ModelSeven/Model.h
        impl/model/ModelSeven/config.h
        impl/model/ModelSeven/SampleScheduler.h
        impl/model/ModelSeven/StateLogger.cpp
        impl/model/ModelSeven/StateLogger.h
        impl/model/ModelSeven/ModelLogger.cpp
        impl/model/ModelSeven/ModelLogger.h
        impl/model/ModelSeven/SamplerLogger.cpp
        impl/model/ModelSeven/SamplerLogger.h

        impl/model/ModelEight/State.cpp
        impl/model/ModelEight/State.h
        impl/model/ModelEight/Model.cpp
        impl/model/ModelEight/Model.h
        impl/model/ModelEight/config.h
        impl/model/ModelEight/SampleScheduler.h
        impl/model/ModelEight/StateLogger.cpp
        impl/model/ModelEight/StateLogger.h
        impl/model/ModelEight/ModelLogger.cpp
        impl/model/ModelEight/ModelLogger.h
        impl/model/ModelEight/SamplerLogger.cpp
        impl/model/ModelEight/SamplerLogger.h

        model/transmission_process/node_transmission_process/SuperInfectionNoMutation.h
        model/observation_process/ObservationProcessLikelihoodv2.h model/transmission_process/node_transmission_process/SimpleLossMutation.h model/transmission_process/node_transmission_process/SimpleLoss.h)


add_library(transmission_networks SHARED STATIC ${SOURCE_FILES})
target_link_libraries(transmission_networks PRIVATE nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX Eigen3::Eigen fmt::fmt mpfr)


install(TARGETS transmission_networks DESTINATION ${TRANSMISSION_NETWORKS_INSTALL_BIN_DIR})