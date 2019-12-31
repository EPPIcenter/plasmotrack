cmake_minimum_required(VERSION 3.15)
project(transmission_networks)

add_subdirectory(core)

set(SOURCE_FILES main.cpp)

include_directories(${PROJECT_BINARY_DIR})

add_executable(transmission_networks ${SOURCE_FILES})
target_link_libraries(transmission_networks transmission_network_core)
install(TARGETS transmission_networks DESTINATION ${TRANSMISSION_NETWORKS_INSTALL_BIN_DIR})