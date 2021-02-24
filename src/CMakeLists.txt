cmake_minimum_required(VERSION 3.17)
project(transmission_networks)


find_package(Boost 1.74 REQUIRED container math_c99l random filesystem)
find_package(Eigen3 REQUIRED)
find_package(nlohmann_json 3.2.0 REQUIRED)
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()



include_directories(SYSTEM ${Boost_INCLUDE_DIR})

set(EIGEN3_INCLUDE_DIR $ENV{EIGEN3_INCLUDE_DIR})
message("Loading Eigen from $ENV{EIGEN3_INCLUDE_DIR}")
if (NOT EIGEN3_INCLUDE_DIR)
    message(FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
endif ()
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

set(SOURCE_FILES
        core/computation/PartialLikelihood.cpp
        core/computation/LogLikelihood.cpp
        core/containers/Locus.cpp
        core/utils/CombinationIndicesGenerator.cpp
        model/observation_process/AlleleCounts.cpp
        core/samplers/general/ConstrainedDiscreteRandomWalk.h
        core/samplers/genetics/RandomAllelesBitSetSampler.h
        core/utils/RandomSequence.h
        core/datatypes/Simplex.cpp
        core/samplers/Scheduler.cpp
        core/samplers/Scheduler.h
        impl/state/ModelOneState.cpp
        impl/state/ModelOneState.h
        impl/model/ModelOne.cpp
        impl/model/ModelOne.h
        core/computation/Transformers/Transformer.h
        core/computation/Transformers/LogTransformer.h
        core/utils/CombinationsWithRepetitionsGenerator.cpp
        core/utils/CombinationsWithRepetitionsGenerator.h
        core/priors/Prior.h
        core/io/parse_json.h
        core/io/loggers/FileOutput.h
        core/io/Loggers/ValueLogger.h
        core/io/serialize.h
        core/io/loggers/FileOutput.cpp
        core/samplers/RandomizedScheduler.h
        core/io/path_parsing.cpp
        core/io/path_parsing.h
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
        core/io/Loggers/LambdaLogger.h
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
        core/io/loggers/AbstractLogger.h core/io/loggers/AbstractOutput.h)


add_library(transmission_networks SHARED STATIC ${SOURCE_FILES})
target_link_libraries(transmission_networks PRIVATE nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)
#target_link_libraries(transmission_networks PRIVATE )
#target_link_libraries(transmission_networks)

#target_link_libraries(transmission_networks transmission_network_model)

install(TARGETS transmission_networks DESTINATION ${TRANSMISSION_NETWORKS_INSTALL_BIN_DIR})