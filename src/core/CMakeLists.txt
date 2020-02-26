cmake_minimum_required(VERSION 3.15)
project(transmission_network_core)

find_package(Boost 1.72 REQUIRED container math_c99l)
find_package(Eigen3 REQUIRED)

set(SOURCE_FILES
        main.cpp
        computation/PartialLikelihood.cpp
        computation/LogLikelihood.cpp
        computation/observation_process/AlleleCounts.cpp
        containers/Locus.cpp
        )

add_library(transmission_network_core SHARED STATIC ${SOURCE_FILES})
set_target_properties(transmission_network_core PROPERTIES LINKER_LANGUAGE CXX)

install(TARGETS transmission_network_core DESTINATION ${TRANSMISSION_NETWORKS_INSTALL_LIB_DIR})