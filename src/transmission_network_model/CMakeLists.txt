cmake_minimum_required(VERSION 3.15)
project(transmission_network_model C CXX)

if(NOT DEFINED transmission_network_model_MAX_ALLELES)
    set(transmission_network_model_MAX_ALLELES 32)
endif()

if(NOT DEFINED transmission_network_model_MAX_GENERATIONS)
    set(transmission_network_model_MAX_GENERATIONS 10)
endif()

if(NOT DEFINED transmission_network_model_MAX_COI)
    set(transmission_network_model_MAX_COI 20)
endif()

configure_file (
    "${PROJECT_SOURCE_DIR}/TransmissionNetworkModelConfig.h.in"
    "${PROJECT_SOURCE_DIR}/TransmissionNetworkModelConfig.h"
)

# Declare Dependencies
find_package(Eigen3 3.3 REQUIRED)
find_package(Boost 1.71 REQUIRED COMPONENTS math)


set(SOURCE_FILES
        LikelihoodFunction.h
        TransmissionMechanism.cpp TransmissionMechanism.h
        TransmissionMechanismParameters.cpp TransmissionMechanismParameters.h
        ObservationMechanism.cpp ObservationMechanism.h
        ObservationMechanismParameters.cpp ObservationMechanismParameters.h
        ModelParameter.cpp ModelParameter.h)

add_library(transmission_network_model SHARED STATIC ${SOURCE_FILES})

include_directories(${PROJECT_BINARY_DIR})

install(TARGETS transmission_network_model DESTINATION ${TRANSMISSION_NETWORKS_INSTALL_LIB_DIR})
install(FILES TransmissionNetworkModel.h DESTINATION ${TRANSMISSION_NETWORKS_INSTALL_INCLUDE_DIR})