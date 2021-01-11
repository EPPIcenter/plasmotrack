cmake_minimum_required(VERSION 3.15)
project(transmission_network_core)

find_package(Boost 1.74 REQUIRED container math_c99l random filesystem)
find_package(Eigen3 REQUIRED)
find_package(GSL REQUIRED)
find_package(nlohmann_json 3.2.0 REQUIRED)
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


set(SOURCE_FILES
        computation/PartialLikelihood.cpp
        computation/LogLikelihood.cpp
        containers/Locus.cpp
        )

add_library(transmission_network_core SHARED STATIC ${SOURCE_FILES})
set_target_properties(transmission_network_core PROPERTIES LINKER_LANGUAGE CXX)

install(TARGETS transmission_network_core DESTINATION ${TRANSMISSION_NETWORKS_INSTALL_LIB_DIR})