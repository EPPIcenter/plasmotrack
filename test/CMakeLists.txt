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


add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
        ${CMAKE_BINARY_DIR}/googletest-build
        EXCLUDE_FROM_ALL)

include_directories(${TRANSMISSION_NETWORK_HEADERS_DIR})

set(SOURCE_FILES
        main.cpp
        src/core/core_likelihood_tests.cpp
        src/core/datatypes/AllelesTest.cpp
        src/core/parameters/OrderingTest.cpp
        src/core/containers/InfectionTest.cpp
        src/core/containers/AlleleFrequencyContainerTest.cpp
        src/core/utils/CombinationIndicesGeneratorTests.cpp
        src/core/datatypes/AlleleFrequenciesVectorTest.cpp
        src/model/OrderDerivedParentSetTest.cpp
        src/model/transmission_process/source_transmission_process/GeometricCOIProbabilityTest.cpp
        src/model/transmission_process/source_transmission_process/MultinomialSourceTransmissionProcessTest.cpp
        src/model/transmission_process/OrderBasedTransmissionProcessTest.cpp)

add_executable(transmission_networks_tests ${SOURCE_FILES})
target_link_libraries(transmission_networks_tests transmission_networks gtest)
install(TARGETS transmission_networks_tests DESTINATION ${TRANSMISSION_NETWORKS_INSTALL_BIN_DIR})