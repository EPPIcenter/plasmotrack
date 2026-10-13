add_executable(direction_check direction_check.cpp)
target_include_directories(direction_check PRIVATE ${TRANSMISSION_NETWORK_HEADERS_DIR})
target_link_libraries(direction_check PRIVATE transmission_networks Boost::boost Eigen3::Eigen nlohmann_json::nlohmann_json fmt::fmt ZLIB::ZLIB)
link_openmp_if_available(direction_check SUPPRESS_PRAGMA_WARNINGS)


add_executable(transplant transplant.cpp)
target_include_directories(transplant PRIVATE ${TRANSMISSION_NETWORK_HEADERS_DIR})
target_link_libraries(transplant PRIVATE transmission_networks Boost::boost Eigen3::Eigen nlohmann_json::nlohmann_json fmt::fmt ZLIB::ZLIB)
link_openmp_if_available(transplant SUPPRESS_PRAGMA_WARNINGS)

add_executable(transplant_groups transplant_groups.cpp)
target_include_directories(transplant_groups PRIVATE ${TRANSMISSION_NETWORK_HEADERS_DIR})
target_link_libraries(transplant_groups PRIVATE transmission_networks Boost::boost Eigen3::Eigen nlohmann_json::nlohmann_json fmt::fmt ZLIB::ZLIB)
link_openmp_if_available(transplant_groups SUPPRESS_PRAGMA_WARNINGS)
