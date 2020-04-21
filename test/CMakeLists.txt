cmake_minimum_required(VERSION 3.15)
project(transmission_networks_tests)

configure_file(${CMAKE_SOURCE_DIR}/cmake/GoogleTest-CMakeLists.txt.in ${CMAKE_BINARY_DIR}/googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

find_package(Boost 1.72 REQUIRED container math_c99l random)
find_package(Eigen3 REQUIRED)

include_directories(SYSTEM ${Boost_INCLUDE_DIR})

set(EIGEN3_INCLUDE_DIR $ENV{EIGEN3_INCLUDE_DIR})
if( NOT EIGEN3_INCLUDE_DIR )
    message( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
endif()
include_directories (SYSTEM ${EIGEN3_INCLUDE_DIR})

add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
        ${CMAKE_BINARY_DIR}/googletest-build
        EXCLUDE_FROM_ALL)

include_directories(${TRANSMISSION_NETWORK_HEADERS_DIR})

set(SOURCE_FILES
        main.cpp
        src/core/core_likelihood_tests.cpp
        src/core/containers/InfectionTest.cpp
        src/core/containers/AlleleFrequencyContainerTest.cpp
        src/core/datatypes/AllelesTest.cpp
        src/core/datatypes/SimplexTest.cpp
        src/core/parameters/OrderingTest.cpp
        src/core/samplers/ContinuousRandomWalkTest.cpp
        src/core/samplers/ConstrainedContinuousRandomWalkTest.cpp
        src/core/samplers/SALTSamplerTest.cpp
        src/core/samplers/OrderSamplerTest.cpp
        src/core/samplers/DiscreteRandomWalkTest.cpp
        src/core/samplers/ConstrainedDiscreteRandomWalkTest.cpp
        src/core/utils/CombinationIndicesGeneratorTests.cpp
        src/model/transmission_process/OrderDerivedParentSetTest.cpp
        src/model/transmission_process/OrderBasedTransmissionProcessTest.cpp
        src/model/transmission_process/source_transmission_process/GeometricCOIProbabilityTest.cpp
        src/model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcessTest.cpp
        src/model/transmission_process/node_transmission_process/GeometricGenerationProbabilityTest.cpp
        src/core/samplers/genetics/RandomAllelesBitSetSamplerTest.cpp
        src/core/utils/RandomSequenceTests.cpp src/model/observation_process/ObservationProcessTest.cpp src/impl/ModelOneTest.cpp src/model/transmission_process/node_transmission_process/ZTMultiplicativeBinomialTest.cpp)

add_executable(transmission_networks_tests ${SOURCE_FILES})
target_link_libraries(transmission_networks_tests transmission_networks gtest)
install(TARGETS transmission_networks_tests DESTINATION ${TRANSMISSION_NETWORKS_INSTALL_BIN_DIR})