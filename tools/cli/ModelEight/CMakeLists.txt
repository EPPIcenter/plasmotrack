cmake_minimum_required(VERSION 3.17)
include(FetchContent)

project(transmission_networks_model_eight)

#FetchContent_Declare(json
#        GIT_REPOSITORY https://github.com/nlohmann/json.git
#        GIT_TAG master
#        GIT_SHALLOW TRUE
#        GIT_PROGRESS TRUE)
#FetchContent_MakeAvailable(json)
#
#FetchContent_Declare(
#        Eigen
#        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
#        GIT_TAG master
#        GIT_SHALLOW TRUE
#        GIT_PROGRESS TRUE)
#set(EIGEN_BUILD_DOC OFF)
#set(BUILD_TESTING OFF)
#set(EIGEN_BUILD_PKGCONFIG OFF)
#FetchContent_MakeAvailable(Eigen)


find_package(Boost)
find_package(Eigen3)
find_package(nlohmann_json)
find_package(OpenMP)

if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

#include_directories(SYSTEM ${Boost_INCLUDE_DIR})

include_directories(${TRANSMISSION_NETWORK_HEADERS_DIR})

set(SOURCE_FILES ModelEight.cpp)
add_executable(transmission_networks_model_eight ${SOURCE_FILES})
target_link_libraries(transmission_networks_model_eight

        PRIVATE
        transmission_networks
        nlohmann_json::nlohmann_json
        Eigen3::Eigen
        OpenMP::OpenMP_CXX
        ${Boost_LIBRARIES}
        )

message(STATUS ${TRANSMISSION_NETWORKS_INSTALL_BIN_DIR})
install(TARGETS transmission_networks_model_eight DESTINATION ${TRANSMISSION_NETWORKS_INSTALL_BIN_DIR})