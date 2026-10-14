include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        main
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

set(TEST_MAIN
    main.cpp
)
set(CORE_CONTAINERS_TESTS
    src/core/containers/InfectionTest.cpp
    src/core/containers/UidOrderTest.cpp
    src/core/containers/AlleleFrequencyContainerTest.cpp
)

set(CORE_DATATYPES_TESTS
    src/core/datatypes/AllelesTest.cpp
    src/core/datatypes/SimplexTest.cpp
)

set(CORE_DISTRIBUTIONS_TESTS
    src/core/distributions/GeometricTest.cpp
    src/core/distributions/DiscreteTest.cpp
    src/core/distributions/ZTGeometricTest.cpp
)

set(CORE_SAMPLERS_TESTS
    src/core/samplers/ContinuousRandomWalkTest.cpp
    src/core/samplers/ConstrainedContinuousRandomWalkTest.cpp
    src/core/samplers/SALTSamplerTest.cpp
    src/core/samplers/genetics/RandomAllelesBitSetSampler4Test.cpp
    src/core/samplers/genetics/GenotypeSamplerStationaryTest.cpp
    src/core/samplers/specialized/CloneOrderSwapSamplerTest.cpp
    src/core/samplers/specialized/InformedOrderJumpSamplerTest.cpp
    src/core/samplers/genetics/StrainCountSamplerTest.cpp
    src/core/samplers/SchedulerTest.cpp
    src/core/samplers/ProposalTest.cpp
)

set(CORE_UTILS_TESTS
    src/core/utils/CombinationIndicesGeneratorTests.cpp
    src/core/utils/InterruptWatcherTest.cpp
    src/core/utils/RandomSequenceTests.cpp
    src/core/utils/CombinationsWithRepetitionsGeneratorTest.cpp
    src/core/utils/ParseJSONTests.cpp
    src/core/utils/ProbAnyMissingTests.cpp
    src/core/utils/NumericsTest.cpp
)

set(CORE_COMPUTATION_TESTS
    src/core/computation/ObservationTimeDerivedOrderingTest.cpp
    src/core/parameters/OrderingTest.cpp
)

set(MODEL_TRANSMISSION_TESTS
    src/model/transmission_process/OrderDerivedParentSetTest.cpp
    src/model/transmission_process/node_transmission_process/MultinomialTransmissionProcessTest.cpp
    src/model/transmission_process/node_transmission_process/MarginalLatentParentTest.cpp
    src/impl/ModelCacheConsistencyTest.cpp
    src/impl/RunRecorderTest.cpp
    src/impl/RegisteredSamplerTest.cpp
    src/impl/QuantityTest.cpp
    src/impl/ReplicaExchangeTest.cpp
    src/impl/CheckpointTest.cpp
    src/impl/RunDriverTest.cpp
    src/impl/ConsoleProgressTest.cpp
    src/impl/RunRecipeTest.cpp
    src/core/samplers/SamplingProbTest.cpp
)

set(MODEL_OBSERVATION_TESTS
    src/model/observation_process/ObservationProcessTest.cpp
)

set(TOOLS_TESTS
    src/tools/GoldenRunTest.cpp
    src/tools/CliStoreTest.cpp
    src/tools/CliCheckpointTest.cpp
    src/tools/CliFailureTest.cpp
    src/core/io/ParquetRoundTripTest.cpp
    src/core/io/ParquetStoreSinkTest.cpp
    src/core/io/DurableFileTest.cpp
)

set(SOURCE_FILES
    ${TEST_MAIN}
    ${CORE_CONTAINERS_TESTS}
    ${CORE_DATATYPES_TESTS}
    ${CORE_DISTRIBUTIONS_TESTS}
    ${CORE_SAMPLERS_TESTS}
    ${CORE_UTILS_TESTS}
    ${CORE_COMPUTATION_TESTS}
    ${MODEL_TRANSMISSION_TESTS}
    ${MODEL_OBSERVATION_TESTS}
    ${TOOLS_TESTS}
)

add_executable(transmission_networks_tests ${SOURCE_FILES})

set_target_properties(transmission_networks_tests PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)
target_compile_features(transmission_networks_tests PRIVATE cxx_std_20)

target_include_directories(transmission_networks_tests
    PRIVATE
        ${TRANSMISSION_NETWORK_HEADERS_DIR}
)

target_link_libraries(transmission_networks_tests
    PRIVATE
        transmission_networks
        Boost::boost
        nlohmann_json::nlohmann_json
        Eigen3::Eigen
        fmt::fmt
        ZLIB::ZLIB
        Parquet::parquet_static
        gtest
)

# The golden-run check runs the CLI, so build it first and tell the test where it and its data are.
add_dependencies(transmission_networks_tests transmission_networks_model)
target_compile_definitions(transmission_networks_tests
    PRIVATE
        GOLDEN_CLI_PATH="$<TARGET_FILE:transmission_networks_model>"
        GOLDEN_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/golden"
        EXAMPLE_INPUT_PATH="${PROJECT_SOURCE_DIR}/examples/example_input.json"
)

add_test(NAME transmission_networks_tests COMMAND transmission_networks_tests)

link_openmp_if_available(transmission_networks_tests)

include(CheckIPOSupported)
check_ipo_supported(RESULT supported OUTPUT error)

if (supported)
    message(STATUS "IPO / LTO enabled")
    set_property(TARGET transmission_networks_tests PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else ()
    message(STATUS "IPO / LTO not supported: <${error}>")
endif ()

install(TARGETS transmission_networks_tests
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    COMPONENT tests
)
