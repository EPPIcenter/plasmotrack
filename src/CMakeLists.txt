cmake_minimum_required(VERSION 3.15)
project(transmission_networks)

find_package(Boost 1.72 REQUIRED container math_c99l)
find_package(Eigen3 REQUIRED)
find_package(GSL REQUIRED)


set(SOURCE_FILES
        core/computation/PartialLikelihood.cpp
        core/computation/LogLikelihood.cpp
        core/containers/Locus.cpp
        core/utils/CombinationIndicesGenerator.cpp
        model/observation_process/AlleleCounts.cpp
        )


add_library(transmission_networks SHARED STATIC ${SOURCE_FILES})
#target_link_libraries(transmission_networks transmission_network_model)

install(TARGETS transmission_networks DESTINATION ${TRANSMISSION_NETWORKS_INSTALL_BIN_DIR})