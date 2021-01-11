cmake_minimum_required(VERSION 3.17)
project(transmission_networks_app)

# App Compile Settings
set(transmission_networks_VERSION_MAJOR 1)
set(transmission_networks_VERSION_MINOR 0)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -pedantic-errors")
if ( CMAKE_COMPILER_IS_GNUCC )
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g -fno-omit-frame-pointer")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -g -march=native -frename-registers -funroll-loops -fno-omit-frame-pointer") # -Ofast is throwing nan errors
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -pedantic-errors")
#set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native -pedantic-errors")
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message("Compiling Clang")
    set(Clang_DIR "/usr/local/opt/llvm/lib/cmake/clang")
    find_package(Clang)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g -fno-omit-frame-pointer")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Ofast -g -march=native -ffast-math")
endif()

#if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
#    message("Compiling Intel")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++ -std=c++2a -gcc-name=gcc-10 -gxx-name=g++-10 -no-use-clang-env")
#    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g")
#    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Ofast -g -march=native -ffast-math")
##    set(CMAKE_CXX)
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