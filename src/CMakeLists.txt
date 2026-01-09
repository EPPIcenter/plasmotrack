set(SPLINE_PATH ${TRANSMISSION_NETWORKS_LIB_DIR}/spline)
message(STATUS "Spline path: ${SPLINE_PATH}")

set(SPLINE_HEADERS
    ${SPLINE_PATH}/spline.h
)

add_library("spline" INTERFACE ${SPLINE_HEADERS})
target_include_directories("spline" INTERFACE
    $<BUILD_INTERFACE:${SPLINE_PATH}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include/spline>
)

set_target_properties("spline" PROPERTIES
    EXPORT_NAME spline
)

set(CORE_COMPUTATION_SOURCES
    core/computation/PartialLikelihood.cpp
    core/computation/LogLikelihood.cpp
    core/computation/ConstrainedOrderDerivedParentSet.h
    core/computation/ObservationTimeDerivedOrdering.h
    core/computation/transformers/Transformer.h
    core/computation/transformers/LogTransformer.h
    core/computation/transformers/Tempered.h
)

set(CORE_CONTAINERS_SOURCES
    core/containers/Locus.cpp
)

set(CORE_DATATYPES_SOURCES
    core/datatypes/Simplex.cpp
    core/datatypes/AllelesBitVec.h
)

set(CORE_DISTRIBUTIONS_SOURCES
    core/distributions/ZTPoisson.h
    core/distributions/DiscreteDistribution.h
    core/distributions/DiscreteDistribution.cpp
    core/distributions/pdfs/DiscretePDF.h
    core/distributions/pdfs/DiscretePDF.cpp
    core/distributions/pdfs/BetaLogPDF.h
    core/distributions/pdfs/BetaLogPDF.cpp
    core/distributions/pdfs/GammaLogPDF.cpp
    core/distributions/pdfs/GammaLogPDF.h
)

set(CORE_IO_SOURCES
    core/io/parse_json.h
    core/io/serialize.h
    core/io/utils.cpp
    core/io/utils.h
    core/io/loggers/FileOutput.h
    core/io/loggers/FileOutput.cpp
    core/io/loggers/ValueLogger.h
    core/io/loggers/MultiValueLogger.h
    core/io/loggers/LambdaLogger.h
    core/io/loggers/ParentSetDistLogger.h
    core/io/loggers/AbstractLogger.h
    core/io/loggers/AbstractOutput.h
    core/io/loggers/CompressedFileOutput.cpp
    core/io/loggers/CompressedFileOutput.h
)

set(CORE_MCMC_SOURCES
    core/mcmc/Chain.h
)

set(CORE_PRIORS_SOURCES
    core/priors/Prior.h
)

set(CORE_SAMPLERS_SOURCES
    core/samplers/SamplerEnum.h
    core/samplers/general/ConstrainedDiscreteRandomWalk.h
    core/samplers/general/SimplexSampler.h
    core/samplers/genetics/RandomAllelesBitSetSampler.h
    core/samplers/genetics/RandomAllelesBitSetSampler4.h
    core/samplers/genetics/ZanellaAllelesBitSetSampler.h
    core/samplers/scheduler/Scheduler.cpp
    core/samplers/scheduler/Scheduler.h
    core/samplers/scheduler/RandomizedScheduler.h
    core/samplers/specialized/ZanellaJointGeneticsGraphSampler.h
    core/samplers/specialized/JointGeneticsTimeSampler.h
    core/samplers/topology/RandomAddEdgeSampler.h
    core/samplers/topology/RandomRemoveEdgeSampler.h
    core/samplers/topology/RandomSwapEdgeSampler.h
    core/samplers/topology/RandomReverseEdgeSampler.h
    core/samplers/topology/ZanellaNeighborOrderSampler.h
    core/samplers/topology/ZanellaOrderSampler.h
    core/samplers/meta/ReplicaExchange.h
)

set(CORE_UTILS_SOURCES
    core/utils/generators/CombinationIndicesGenerator.cpp
    core/utils/generators/CombinationsWithRepetitionsGenerator.cpp
    core/utils/generators/CombinationsWithRepetitionsGenerator.h
    core/utils/generators/RandomSequence.h
    core/utils/ProbAnyMissing.h
    core/utils/ProbAnyMissing.cpp
    core/utils/LogPQ.cpp
    core/utils/LogPQ.h
    core/utils/timers.h
)

set(MODEL_SOURCES
    model/OrderingLikelihood.h
    model/observation_process/AlleleCounts.cpp
    model/observation_process/ObservationProcessLikelihoodv2.h
    model/transmission_process/NetworkBasedTransmissionProcess.h
    model/transmission_process/node_transmission_process/NoSuperInfectionMutation.h
    model/transmission_process/node_transmission_process/SuperInfectionNoMutation.h
    model/transmission_process/node_transmission_process/SimpleLossMutation.h
    model/transmission_process/node_transmission_process/SimpleLoss.h
    model/transmission_process/node_transmission_process/MultinomialTransmissionProcess.h
    model/transmission_process/node_transmission_process/MultinomialTransmissionProcess2.h
    model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcess.h
)

set(IMPL_SOURCES
    impl/model/Model/State.cpp
    impl/model/Model/State.h
    impl/model/Model/Model.cpp
    impl/model/Model/Model.h
    impl/model/Model/config.h
    impl/model/Model/ModelLogger.cpp
    impl/model/Model/ModelLogger.h
    impl/model/Model/SampleScheduler.h
    impl/model/Model/SequentialScheduler.h
    impl/model/Model/StateLogger.cpp
    impl/model/Model/StateLogger.h
)

set(OTHER_SOURCES
    core/config.h
    core/containers/AllowedRelationships.h
    ../lib/spline/spline.h
)

set(SOURCE_FILES
    ${CORE_COMPUTATION_SOURCES}
    ${CORE_CONTAINERS_SOURCES}
    ${CORE_DATATYPES_SOURCES}
    ${CORE_DISTRIBUTIONS_SOURCES}
    ${CORE_IO_SOURCES}
    ${CORE_MCMC_SOURCES}
    ${CORE_PRIORS_SOURCES}
    ${CORE_SAMPLERS_SOURCES}
    ${CORE_UTILS_SOURCES}
    ${MODEL_SOURCES}
    ${IMPL_SOURCES}
    ${OTHER_SOURCES}
)

add_library(transmission_networks STATIC ${SOURCE_FILES})

set_target_properties(transmission_networks PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)
target_compile_features(transmission_networks PUBLIC cxx_std_20)

target_include_directories(transmission_networks
    PUBLIC
        $<BUILD_INTERFACE:${TRANSMISSION_NETWORK_HEADERS_DIR}>
        $<BUILD_INTERFACE:${TRANSMISSION_NETWORKS_INSTALL_INCLUDE_DIR}>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
    PRIVATE
        ${TRANSMISSION_NETWORKS_LIB_DIR}
)

target_link_libraries(transmission_networks
    PUBLIC
        spline
    PRIVATE
        Boost::boost
        Eigen3::Eigen
        nlohmann_json::nlohmann_json
        fmt::fmt
        ZLIB::ZLIB
)

link_openmp_if_available(transmission_networks)

install(TARGETS transmission_networks spline
    EXPORT transmission_networksTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        COMPONENT transmission_networks
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        COMPONENT transmission_networks
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
        COMPONENT transmission_networks
    INCLUDES DESTINATION ${CMAKE_INSTALL_PREFIX}/include
        COMPONENT transmission_networks
)

install(DIRECTORY ${TRANSMISSION_NETWORK_HEADERS_DIR}/
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include
    COMPONENT transmission_networks
    EXCLUDE_FROM_ALL
    FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
)

install(FILES ${SPLINE_HEADERS}
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/spline
    COMPONENT transmission_networks
)

install(EXPORT transmission_networksTargets
    FILE transmission_networksTargets.cmake
    NAMESPACE transmission_networks::
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/transmission_networks
    COMPONENT transmission_networks
    EXPORT_LINK_INTERFACE_LIBRARIES
)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/transmission_networksConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
    ${CMAKE_SOURCE_DIR}/cmake/transmission_networksConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/transmission_networksConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/transmission_networks
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/transmission_networksConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/transmission_networksConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/transmission_networks
    COMPONENT transmission_networks
)
