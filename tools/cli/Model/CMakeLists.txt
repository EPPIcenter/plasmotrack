set(SOURCE_FILES Model.cpp)

add_executable(transmission_networks_model ${SOURCE_FILES})

set_target_properties(transmission_networks_model PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)
target_compile_features(transmission_networks_model PRIVATE cxx_std_20)

target_include_directories(transmission_networks_model
    PRIVATE
        ${TRANSMISSION_NETWORK_HEADERS_DIR}
)

target_link_libraries(transmission_networks_model
    PRIVATE
        transmission_networks
        Boost::program_options
        Eigen3::Eigen
        nlohmann_json::nlohmann_json
        fmt::fmt
        ZLIB::ZLIB
)

# Configure version information
# Ensure PATCH version is always numeric (default to 0 if empty)
if("${PROJECT_VERSION_PATCH}" STREQUAL "")
    set(VERSION_PATCH 0)
else()
    set(VERSION_PATCH ${PROJECT_VERSION_PATCH})
endif()

target_compile_definitions(transmission_networks_model
    PRIVATE
        TRANSMISSION_NETWORKS_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
        TRANSMISSION_NETWORKS_VERSION_MINOR=${PROJECT_VERSION_MINOR}
        TRANSMISSION_NETWORKS_VERSION_PATCH=${VERSION_PATCH}
)

link_openmp_if_available(transmission_networks_model SUPPRESS_PRAGMA_WARNINGS)

install(TARGETS transmission_networks_model
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    COMPONENT cli
)
