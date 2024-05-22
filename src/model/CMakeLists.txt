cmake_minimum_required(VERSION 3.15)
project(transmission_network_model)


set(
        SOURCE_FILES
        observation_process/AlleleCounts.cpp
)

add_library(transmission_network_model SHARED STATIC ${SOURCE_FILES})
set_target_properties(transmission_network_model PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(transmission_network_model transmission_network_core)

install(TARGETS transmission_network_model DESTINATION ${TRANSMISSION_NETWORKS_INSTALL_LIB_DIR})