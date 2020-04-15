cmake_minimum_required(VERSION 3.15)
project(transmission_networks)


find_package(Boost 1.72 REQUIRED container math_c99l random)
find_package(Eigen3 REQUIRED)

include_directories(SYSTEM ${Boost_INCLUDE_DIR})

set(EIGEN3_INCLUDE_DIR $ENV{EIGEN3_INCLUDE_DIR})
if( NOT EIGEN3_INCLUDE_DIR )
    message( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
endif()
include_directories (${EIGEN3_INCLUDE_DIR})

set(SOURCE_FILES
        core/computation/PartialLikelihood.cpp
        core/computation/LogLikelihood.cpp
        core/containers/Locus.cpp
        core/utils/CombinationIndicesGenerator.cpp
        model/observation_process/AlleleCounts.cpp
        core/samplers/OrderSampler.h)


add_library(transmission_networks SHARED STATIC ${SOURCE_FILES})
target_link_libraries(transmission_networks)

#target_link_libraries(transmission_networks transmission_network_model)

install(TARGETS transmission_networks DESTINATION ${TRANSMISSION_NETWORKS_INSTALL_BIN_DIR})