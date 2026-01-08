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

link_openmp_if_available(transmission_networks_model SUPPRESS_PRAGMA_WARNINGS)

install(TARGETS transmission_networks_model
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    COMPONENT cli
)
