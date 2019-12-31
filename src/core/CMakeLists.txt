cmake_minimum_required(VERSION 3.15)
project(transmission_network_core)

set(SOURCE_FILES
        Likelihood/Likelihood.cpp
        Likelihood/PartialLikelihood.cpp
        Likelihood/AdderLikelihood.cpp
        )

add_library(transmission_network_core SHARED STATIC ${SOURCE_FILES})
set_target_properties(transmission_network_core PROPERTIES LINKER_LANGUAGE CXX)

include_directories(${PROJECT_BINARY_DIR})

install(TARGETS transmission_network_core DESTINATION ${TRANSMISSION_NETWORKS_INSTALL_LIB_DIR})