cmake_minimum_required(VERSION 3.15)
project(transmission_networks_app)

# App Compile Settings
set(transmission_networks_VERSION_MAJOR 1)
set(transmission_networks_VERSION_MINOR 0)

set(CMAKE_CXX_STANDARD 17)


if ( CMAKE_COMPILER_IS_GNUCC )
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g -Werror -pedantic-errors")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Ofast -g -march=native -ffast-math -Werror -pedantic-errors -fno-omit-frame-pointer")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -pedantic-errors")
#set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native -pedantic-errors")
endif()

#if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#    message("Compiling Clang")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mlinker-version=305 -Wall -Wextra")
#    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g -Werror -pedantic-errors")
#    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Ofast -march=native -ffast-math -Werror -pedantic-errors")
#endif()

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

set(TRANSMISSION_NETWORKS_INSTALL_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(TRANSMISSION_NETWORKS_INSTALL_BIN_DIR ${PROJECT_SOURCE_DIR}/bin)
set(TRANSMISSION_NETWORKS_LIB_DIR ${PROJECT_SOURCE_DIR}/lib)

set(TRANSMISSION_NETWORK_HEADERS_DIR ${PROJECT_SOURCE_DIR}/src)

include_directories(${TRANSMISSION_NETWORKS_INSTALL_INCLUDE_DIR})
include_directories(${TRANSMISSION_NETWORK_MODEL_HEADERS_DIR})
include_directories(${TRANSMISSION_NETWORK_HEADERS_DIR})
include_directories(${PROJECT_BINARY_DIR})


add_subdirectory(src)
add_subdirectory(test)