cmake_minimum_required(VERSION 3.17)
project(transmission_networks_model_five)

find_package(Boost 1.74 COMPONENTS program_options REQUIRED container math_c99l random filesystem)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(nlohmann_json 3.2.0 REQUIRED)
find_package(OpenMP)


if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

include_directories(SYSTEM ${Boost_INCLUDE_DIR})

include_directories(${TRANSMISSION_NETWORK_HEADERS_DIR})

set(SOURCE_FILES ModelFive.cpp)
add_executable(transmission_networks_model_five ${SOURCE_FILES})
target_link_libraries(transmission_networks_model_five
        PRIVATE
        transmission_networks
        nlohmann_json::nlohmann_json
        Eigen3::Eigen
        OpenMP::OpenMP_CXX
        ${Boost_LIBRARIES}
        )

message(${TRANSMISSION_NETWORKS_INSTALL_BIN_DIR})
install(TARGETS transmission_networks_model_five DESTINATION ${TRANSMISSION_NETWORKS_INSTALL_BIN_DIR})