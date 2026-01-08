function(find_transmission_networks_dependencies)
    find_package(Boost REQUIRED COMPONENTS program_options)
    find_package(Eigen3 REQUIRED)
    find_package(nlohmann_json REQUIRED)
    find_package(fmt REQUIRED)
    find_package(ZLIB REQUIRED)
    
    find_package(OpenMP COMPONENTS CXX)
    
    if(OpenMP_FOUND AND OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found - parallel computation enabled")
        message(STATUS "  OpenMP CXX flags: ${OpenMP_CXX_FLAGS}")
        message(STATUS "  OpenMP CXX libs: ${OpenMP_CXX_LIB_NAMES}")
        set(OpenMP_FOUND TRUE PARENT_SCOPE)
        set(OpenMP_CXX_FOUND TRUE PARENT_SCOPE)
        if(OpenMP_CXX_FLAGS)
            set(OpenMP_CXX_FLAGS "${OpenMP_CXX_FLAGS}" PARENT_SCOPE)
        endif()
    elseif(OpenMP_FOUND AND NOT OpenMP_CXX_FOUND)
        message(WARNING "OpenMP found but C++ component not available - parallel computation disabled")
        set(OpenMP_FOUND FALSE PARENT_SCOPE)
        set(OpenMP_CXX_FOUND FALSE PARENT_SCOPE)
    else()
        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
            message(STATUS "OpenMP package not found, trying compiler detection...")
            include(CheckCXXCompilerFlag)
            check_cxx_compiler_flag("-fopenmp" COMPILER_SUPPORTS_OPENMP)
            
            if(COMPILER_SUPPORTS_OPENMP)
                message(STATUS "Compiler supports OpenMP, creating manual OpenMP target")
                set(OpenMP_FOUND TRUE PARENT_SCOPE)
                set(OpenMP_CXX_FOUND TRUE PARENT_SCOPE)
                set(OpenMP_CXX_FLAGS "-fopenmp" PARENT_SCOPE)
                set(OpenMP_FOUND TRUE)
                set(OpenMP_CXX_FOUND TRUE)
                set(OpenMP_CXX_FLAGS "-fopenmp")
                
                if(NOT TARGET OpenMP::OpenMP_CXX)
                    add_library(OpenMP::OpenMP_CXX INTERFACE IMPORTED)
                    set_target_properties(OpenMP::OpenMP_CXX PROPERTIES
                        INTERFACE_COMPILE_OPTIONS "-fopenmp"
                        INTERFACE_LINK_OPTIONS "-fopenmp"
                    )
                endif()
                message(STATUS "OpenMP enabled via compiler support (-fopenmp)")
            else()
                message(STATUS "OpenMP not found - parallel computation disabled")
                message(STATUS "  Install libomp-dev (Ubuntu/Debian) or libomp (macOS) to enable OpenMP")
                set(OpenMP_FOUND FALSE PARENT_SCOPE)
                set(OpenMP_CXX_FOUND FALSE PARENT_SCOPE)
            endif()
        else()
            message(STATUS "OpenMP not found - parallel computation disabled")
            message(STATUS "  Install libomp-dev (Ubuntu/Debian) or libomp (macOS) to enable OpenMP")
            set(OpenMP_FOUND FALSE PARENT_SCOPE)
            set(OpenMP_CXX_FOUND FALSE PARENT_SCOPE)
        endif()
    endif()
endfunction()

function(link_openmp_if_available target_name)
    cmake_parse_arguments(ARGS "SUPPRESS_PRAGMA_WARNINGS" "" "" ${ARGN})
    
    if(OpenMP_FOUND AND OpenMP_CXX_FOUND)
        target_link_libraries(${target_name}
            PRIVATE
                OpenMP::OpenMP_CXX
        )
        if(OpenMP_CXX_FLAGS)
            target_compile_options(${target_name}
                PRIVATE
                    $<$<COMPILE_LANGUAGE:CXX>:${OpenMP_CXX_FLAGS}>
            )
            message(STATUS "OpenMP linked to ${target_name} with flags: ${OpenMP_CXX_FLAGS}")
        else()
            if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
                target_compile_options(${target_name}
                    PRIVATE
                        -fopenmp
                )
                message(STATUS "OpenMP linked to ${target_name} with -fopenmp flag")
            endif()
        endif()
    elseif(ARGS_SUPPRESS_PRAGMA_WARNINGS)
        target_compile_options(${target_name}
            PRIVATE
                -Wno-unknown-pragmas
        )
        message(STATUS "OpenMP not available for ${target_name}, suppressing pragma warnings")
    endif()
endfunction()
