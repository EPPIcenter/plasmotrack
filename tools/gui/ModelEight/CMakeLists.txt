cmake_minimum_required(VERSION 3.20)
project(transmission_networks_model_eight_gui)

find_package(Boost)
find_package(Eigen3)
find_package(nlohmann_json)
find_package(OpenMP)

find_package(glfw3)
find_package(OpenGL)


#set(IMGUI_PATH "${TRANSMISSION_NETWORKS_LIB_DIR}/imgui")
#file(GLOB IMGUI_SOURCES ${IMGUI_PATH}/**/*.cpp ${IMGUI_PATH}/*.cpp ${IMGUI_PATH}/*.h ${IMGUI_PATH}/**/*.h)
#add_library("imgui" STATIC ${IMGUI_SOURCES})
#target_include_directories("imgui" PUBLIC ${IMGUI_PATH})

set(IMGUI_PATH ${TRANSMISSION_NETWORKS_LIB_DIR}/imgui)
message(STATUS ${IMGUI_PATH})
file(GLOB IMGUI_GLOB
        ${IMGUI_PATH}/imgui.h
        ${IMGUI_PATH}/imgui.cpp
        ${IMGUI_PATH}/imconfig.h
        ${IMGUI_PATH}/imgui_demo.cpp
        ${IMGUI_PATH}/imgui_draw.cpp
        ${IMGUI_PATH}/imgui_internal.h
        ${IMGUI_PATH}/imstb_rectpack.h
        ${IMGUI_PATH}/imstb_textedit.h
        ${IMGUI_PATH}/imstb_truetype.h
        ${IMGUI_PATH}/imgui_tables.cpp
        ${IMGUI_PATH}/imgui_widgets.cpp

        # specific bindings...
        ${IMGUI_PATH}/backends/imgui_impl_glfw.h
        ${IMGUI_PATH}/backends/imgui_impl_glfw.cpp
        ${IMGUI_PATH}/backends/imgui_impl_opengl3.h
        ${IMGUI_PATH}/backends/imgui_impl_opengl3.cpp
        ${IMGUI_PATH}/backends/imgui_impl_opengl3_loader.cpp)
add_library("imgui" STATIC ${IMGUI_GLOB})
target_include_directories("imgui" PUBLIC ${IMGUI_PATH})
target_link_libraries("imgui" PRIVATE glfw OpenGL::GL)


#find_package(imgui REQUIRED)

if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

include_directories(SYSTEM ${Boost_INCLUDE_DIR})

include_directories(${TRANSMISSION_NETWORK_HEADERS_DIR})

set(SOURCE_FILES main.cpp Application.cpp Application.h)
add_executable(transmission_networks_model_eight_gui ${SOURCE_FILES})
target_link_libraries(transmission_networks_model_eight_gui
        PUBLIC
        imgui
        PRIVATE
        transmission_networks
        nlohmann_json::nlohmann_json
        Eigen3::Eigen
        OpenMP::OpenMP_CXX
        ${Boost_LIBRARIES}
        glfw
        )

message(${TRANSMISSION_NETWORKS_INSTALL_BIN_DIR})
install(TARGETS transmission_networks_model_eight_gui DESTINATION ${TRANSMISSION_NETWORKS_INSTALL_BIN_DIR})
