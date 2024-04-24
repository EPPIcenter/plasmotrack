cmake_minimum_required(VERSION 3.17)

project(transmission_networks_model_nine)

find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(fmt REQUIRED)
find_package(OpenMP)
find_package(ZLIB REQUIRED)

if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

include_directories(SYSTEM ${Boost_INCLUDE_DIR})

include_directories(${TRANSMISSION_NETWORK_HEADERS_DIR})

set(SOURCE_FILES ModelNine.cpp)
add_executable(transmission_networks_model_nine ${SOURCE_FILES})
target_link_libraries(transmission_networks_model_nine
        PRIVATE
        transmission_networks
        nlohmann_json::nlohmann_json
        Eigen3::Eigen
        OpenMP::OpenMP_CXX
        ${Boost_LIBRARIES}
        fmt::fmt
        ZLIB::ZLIB
        )

message(STATUS ${TRANSMISSION_NETWORKS_INSTALL_BIN_DIR})
install(TARGETS transmission_networks_model_nine DESTINATION ${TRANSMISSION_NETWORKS_INSTALL_BIN_DIR})
