cmake_minimum_required(VERSION 4.2.1)

list(APPEND CMAKE_PROJECT_TOP_LEVEL_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/conan_provider.cmake")

project(transmission_networks_app VERSION 1.0 LANGUAGES CXX)

set(CMAKE_FIND_PACKAGE_PREFER_CONFIG TRUE)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_EXPORT_PACKAGE_DEPENDENCIES ON)

option(CMAKE_COLOR_DIAGNOSTICS "Always produce colored output" ON)
option(CMAKE_VERBOSE_MAKEFILE "Enable verbose makefile output" OFF)

if(NOT CMAKE_INSTALL_DEFAULT_COMPONENT_NAME)
    set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME "Unspecified")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
    message(STATUS "Using ccache: ${CCACHE_PROGRAM}")
endif()

option(TRANSMISSION_NETWORKS_WERROR "Treat compiler warnings as errors" OFF)
option(TRANSMISSION_NETWORKS_AGGRESSIVE_OPTIMIZATION "Use -O3 instead of -O2 for Release builds" OFF)
option(TRANSMISSION_NETWORKS_NATIVE_OPTIMIZATION "Use -march=native for Release builds (non-portable)" OFF)

add_compile_options(
    -Wall
    -Wextra
    -Wno-unused-function
    -Wno-unused-parameter
    -Wno-maybe-uninitialized
)

if(TRANSMISSION_NETWORKS_WERROR)
    add_compile_options(-Werror)
    message(STATUS "Warnings as errors enabled")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Debug>>:-Og>
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Debug>>:-fno-omit-frame-pointer>
    )
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    if(TRANSMISSION_NETWORKS_AGGRESSIVE_OPTIMIZATION)
        message(STATUS "Aggressive optimization (-O3) enabled for Release builds")
    else()
        string(REGEX REPLACE "-O[0-9s]" "-O2" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
        if(NOT CMAKE_CXX_FLAGS_RELEASE MATCHES "-O[0-9s]")
            set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
        endif()
        message(STATUS "Using -O2 optimization for Release builds")
    endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options($<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RelWithDebInfo>>:-fno-omit-frame-pointer>)
endif()

if(TRANSMISSION_NETWORKS_NATIVE_OPTIMIZATION)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
        message(STATUS "Native CPU optimization enabled")
    endif()
endif()

include(cmake/Dependencies.cmake)
find_transmission_networks_dependencies()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR} CACHE PATH "Install path prefix" FORCE)
endif()

set(TRANSMISSION_NETWORKS_INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include)
set(TRANSMISSION_NETWORKS_INSTALL_BIN_DIR ${CMAKE_INSTALL_PREFIX}/bin)
set(TRANSMISSION_NETWORKS_LIB_DIR ${PROJECT_SOURCE_DIR}/lib)
set(TRANSMISSION_NETWORK_HEADERS_DIR ${PROJECT_SOURCE_DIR}/src)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(tools/cli/Model)
